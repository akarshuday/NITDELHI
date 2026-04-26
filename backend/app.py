from flask import Flask, request, render_template, jsonify
import os
import sys
import sqlite3
import json
import bcrypt
import jwt
import datetime
import hashlib
import random
import urllib.request
import requests
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from amr_utils import load_model_bundle, parse_fasta, predict_resistance, screen_resistance_markers

# ── Groq API Configuration ──────────────────────────────────────────────────
GROQ_API_KEY = os.getenv("GROQ_API_KEY")
GROQ_BASE_URL = "https://api.groq.com/openai/v1/chat/completions"
GROQ_MODEL = "llama-3.1-8b-instant"

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
MBASE = os.path.join(BASE, '..', 'model')
DB_PATH = os.path.join(BASE, '..', 'results.db')

app = Flask(__name__, template_folder=os.path.join(BASE, '..', 'frontend', 'templates'))
app.config['TEMPLATES_AUTO_RELOAD'] = True

# ── Load Model and Data ──────────────────────────────────────────────────────
bundle = load_model_bundle(os.path.join(BASE, '..'))
feature_cols = bundle['feature_cols']
antibiotic_cols = bundle['antibiotic_cols']
IMPORTANCES = bundle['importances']
GENE_SIGS = bundle['gene_defs']
BASE_FEATURE_COLS = [gene for gene in feature_cols if not gene.endswith('_hits')]

# ── Database Initialization ──────────────────────────────────────────────────
def init_db():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    # Build predictions table dynamically from loaded antibiotic list
    ab_cols = ",\n            ".join(f"{ab.replace('-', '_').replace(' ', '_')} TEXT" for ab in antibiotic_cols)
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS predictions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            filename TEXT,
            {ab_cols},
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    # Add any missing columns for antibiotics added after initial schema creation
    cursor.execute("PRAGMA table_info(predictions)")
    existing_cols = {row[1] for row in cursor.fetchall()}
    for ab in antibiotic_cols:
        col = ab.replace('-', '_').replace(' ', '_')
        if col not in existing_cols:
            cursor.execute(f"ALTER TABLE predictions ADD COLUMN {col} TEXT")
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS users (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            email TEXT UNIQUE NOT NULL,
            password_hash TEXT NOT NULL,
            reset_token TEXT,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    conn.commit()
    conn.close()

init_db()

# ── Load Smart Input Dataset ─────────────────────────────────────────────────
SMART_DATA_PATH = os.path.join(BASE, '..', 'data', 'bacteria_disease_variants.json')
try:
    with open(SMART_DATA_PATH, 'r') as f:
        SMART_DATA = json.load(f)
except Exception as e:
    print(f"Warning: Could not load {SMART_DATA_PATH}: {e}")
    SMART_DATA = {"diseases": {}, "bacteria": {}}

SMART_SEQUENCE_FILLER = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"


def _dedupe_profile_genes(*groups: list[str]) -> list[str]:
    ordered: list[str] = []
    seen = set()
    for group in groups:
        for gene in group:
            if gene in BASE_FEATURE_COLS and gene not in seen:
                ordered.append(gene)
                seen.add(gene)
    return ordered


def _profile_genes(*, exact: tuple[str, ...] = (), prefixes: tuple[str, ...] = (), contains: tuple[str, ...] = ()) -> list[str]:
    genes = list(exact)
    for gene in BASE_FEATURE_COLS:
        if prefixes and any(gene.startswith(prefix) for prefix in prefixes):
            genes.append(gene)
            continue
        if contains and any(token in gene for token in contains):
            genes.append(gene)
    return _dedupe_profile_genes(genes)


AMINOGLYCOSIDE_PROFILE = _profile_genes(
    exact=(
        "16s_rrna_methyltransferase_a1408",
        "16s_rrna_methyltransferase_g1405",
        "aminoglycoside_bifunctional_resistance_protein",
        "capreomycin_phosphotransferase",
    ),
    prefixes=("aac_", "ant_", "aph_"),
)
EFFLUX_AND_PORIN_PROFILE = _profile_genes(
    contains=("efflux_pump", "mate_transporter", "porin", "resistance_nodulation"),
)
COMMON_NON_BETA_PROFILE = _profile_genes(
    exact=(
        "quinolone_resistance_protein_qnr",
        "tetracycline_inactivation_enzyme",
        "tetracycline_resistant_ribosomal_protection_protein",
        "trimethoprim_resistant_dihydrofolate_reductase_dfr",
        "mcr_phosphoethanolamine_transferase",
        "fosfomycin_thiol_transferase",
    ),
)
ENTERIC_BETA_LACTAM_PROFILE = _profile_genes(
    exact=(
        "tem_beta_lactamase",
        "ctx_m_beta_lactamase",
        "shv_beta_lactamase",
        "ampc_type_beta_lactamase",
        "ndm_beta_lactamase",
        "kpc_beta_lactamase",
        "imp_beta_lactamase",
        "vim_beta_lactamase",
        "oxa_beta_lactamase",
    ),
)
NONFERMENTER_BETA_LACTAM_PROFILE = _profile_genes(
    exact=(
        "oxa_beta_lactamase",
        "oxa_beta_lactamase_oxa_23_like_beta_lactamase",
        "oxa_beta_lactamase_oxa_24_like_beta_lactamase",
        "oxa_beta_lactamase_oxa_48_like_beta_lactamase",
        "oxa_beta_lactamase_oxa_50_like_beta_lactamase",
        "oxa_beta_lactamase_oxa_51_like_beta_lactamase",
        "oxa_beta_lactamase_oxa_58_like_beta_lactamase",
        "ndm_beta_lactamase",
        "vim_beta_lactamase",
        "imp_beta_lactamase",
        "ampc_type_beta_lactamase",
    ),
)
ENTERIC_PROFILE = _dedupe_profile_genes(
    ENTERIC_BETA_LACTAM_PROFILE,
    AMINOGLYCOSIDE_PROFILE,
    COMMON_NON_BETA_PROFILE,
    EFFLUX_AND_PORIN_PROFILE,
)
NONFERMENTER_PROFILE = _dedupe_profile_genes(
    NONFERMENTER_BETA_LACTAM_PROFILE,
    AMINOGLYCOSIDE_PROFILE,
    COMMON_NON_BETA_PROFILE,
    EFFLUX_AND_PORIN_PROFILE,
)
STAPH_PROFILE = _dedupe_profile_genes(
    _profile_genes(
        exact=(
            "methicillin_resistant_pbp2",
            "blaz_beta_lactamase",
            "quinolone_resistance_protein_qnr",
            "tetracycline_inactivation_enzyme",
            "tetracycline_resistant_ribosomal_protection_protein",
        ),
    ),
    AMINOGLYCOSIDE_PROFILE,
    _profile_genes(
        exact=("small_multidrug_resistance_smr_antibiotic_efflux_pump",),
    ),
)
ENTEROCOCCUS_PROFILE = _dedupe_profile_genes(
    _profile_genes(
        exact=(
            "aac_6",
            "ant_6",
            "tetracycline_inactivation_enzyme",
            "tetracycline_resistant_ribosomal_protection_protein",
            "trimethoprim_resistant_dihydrofolate_reductase_dfr",
        ),
    ),
    _profile_genes(
        exact=("atp_binding_cassette_abc_antibiotic_efflux_pump",),
    ),
)
RESPIRATORY_PROFILE = _dedupe_profile_genes(
    _profile_genes(
        exact=(
            "tem_beta_lactamase",
            "bro_beta_lactamase",
            "quinolone_resistance_protein_qnr",
            "tetracycline_inactivation_enzyme",
            "tetracycline_resistant_ribosomal_protection_protein",
            "trimethoprim_resistant_dihydrofolate_reductase_dfr",
        ),
    ),
    _profile_genes(
        exact=("general_bacterial_porin_with_reduced_permeability_to_beta_lactams",),
    ),
)
SMART_PROFILE_OVERRIDES = {
    "Acinetobacter baumannii": NONFERMENTER_PROFILE,
    "Enterococcus faecalis": ENTEROCOCCUS_PROFILE,
    "Enterococcus faecium": ENTEROCOCCUS_PROFILE,
    "Escherichia coli": ENTERIC_PROFILE,
    "Klebsiella pneumoniae": ENTERIC_PROFILE,
    "Mycobacterium tuberculosis": ["quinolone_resistance_protein_qnr"],
    "Neisseria gonorrhoeae": ["quinolone_resistance_protein_qnr"],
    "Proteus mirabilis": ENTERIC_PROFILE,
    "Pseudomonas aeruginosa": NONFERMENTER_PROFILE,
    "Salmonella enterica": ENTERIC_PROFILE,
    "Salmonella typhi": ENTERIC_PROFILE,
    "Staphylococcus aureus": STAPH_PROFILE,
    "Staphylococcus saprophyticus": STAPH_PROFILE,
    "Streptococcus pneumoniae": RESPIRATORY_PROFILE,
    "Vibrio cholerae": ENTERIC_PROFILE,
}
SMART_GENUS_PROFILES = {
    "Acinetobacter": NONFERMENTER_PROFILE,
    "Brucella": ["tetracycline_resistant_ribosomal_protection_protein"],
    "Campylobacter": ["quinolone_resistance_protein_qnr", "tetracycline_resistant_ribosomal_protection_protein"],
    "Enterococcus": ENTEROCOCCUS_PROFILE,
    "Escherichia": ENTERIC_PROFILE,
    "Haemophilus": RESPIRATORY_PROFILE,
    "Klebsiella": ENTERIC_PROFILE,
    "Legionella": RESPIRATORY_PROFILE,
    "Moraxella": RESPIRATORY_PROFILE,
    "Mycobacterium": ["quinolone_resistance_protein_qnr"],
    "Mycoplasma": ["quinolone_resistance_protein_qnr"],
    "Neisseria": ["quinolone_resistance_protein_qnr"],
    "Proteus": ENTERIC_PROFILE,
    "Pseudomonas": NONFERMENTER_PROFILE,
    "Rickettsia": ["tetracycline_resistant_ribosomal_protection_protein"],
    "Salmonella": ENTERIC_PROFILE,
    "Shigella": ENTERIC_PROFILE,
    "Staphylococcus": STAPH_PROFILE,
    "Streptococcus": RESPIRATORY_PROFILE,
    "Vibrio": ENTERIC_PROFILE,
}


def _serialize_fasta(header: str, sequence: str, width: int = 80) -> str:
    lines = [header]
    for start in range(0, len(sequence), width):
        lines.append(sequence[start : start + width])
    return "\n".join(lines)


def _default_profile_for_bacteria(bacteria_name: str) -> list[str]:
    digest = hashlib.sha256(bacteria_name.lower().encode("utf-8")).digest()
    genes: list[str] = []
    for index in digest:
        gene = BASE_FEATURE_COLS[index % len(BASE_FEATURE_COLS)]
        if gene not in genes:
            genes.append(gene)
        if len(genes) == 2:
            break
    return genes


def _resolve_profile_for_bacteria(bacteria_name: str) -> list[str]:
    if bacteria_name in SMART_PROFILE_OVERRIDES:
        return [gene for gene in SMART_PROFILE_OVERRIDES[bacteria_name] if gene in BASE_FEATURE_COLS]

    genus = bacteria_name.split(" ", 1)[0]
    if genus in SMART_GENUS_PROFILES:
        return [gene for gene in SMART_GENUS_PROFILES[genus] if gene in BASE_FEATURE_COLS]

    return _default_profile_for_bacteria(bacteria_name)


def _enrich_smart_sequence(bacteria_name: str, raw_sequence: str) -> str:
    sequence = parse_fasta(raw_sequence)
    if not sequence:
        return raw_sequence

    genes = _resolve_profile_for_bacteria(bacteria_name)
    marker_segments = []
    for gene in genes:
        signatures = GENE_SIGS.get(gene, {}).get("signatures", [])
        if signatures:
            marker_segments.append(signatures[0])

    if not marker_segments:
        return raw_sequence

    header = f">{bacteria_name.replace(' ', '_')}_smart_sequence | MicrobeNet"
    enriched_sequence = sequence + SMART_SEQUENCE_FILLER
    for marker in marker_segments:
        if marker not in enriched_sequence:
            enriched_sequence += marker + SMART_SEQUENCE_FILLER
    return _serialize_fasta(header, enriched_sequence)


def _find_bacteria_match(bacteria_name: str) -> tuple[str | None, str | None, list[str]]:
    if not bacteria_name:
        return None, None, []

    bacteria_map = SMART_DATA.get('bacteria', {})
    if bacteria_name in bacteria_map:
        return bacteria_name, bacteria_map[bacteria_name], []

    normalized_query = bacteria_name.strip().lower()
    exact_matches = [name for name in bacteria_map if name.lower() == normalized_query]
    if len(exact_matches) == 1:
        match = exact_matches[0]
        return match, bacteria_map[match], []

    prefix_matches = [name for name in bacteria_map if name.lower().startswith(normalized_query)]
    if len(prefix_matches) == 1:
        match = prefix_matches[0]
        return match, bacteria_map[match], []

    contains_matches = [name for name in bacteria_map if normalized_query in name.lower()]
    if len(contains_matches) == 1:
        match = contains_matches[0]
        return match, bacteria_map[match], []

    candidates = prefix_matches or contains_matches
    return None, None, sorted(candidates)



# ── Authentication Routes ────────────────────────────────────────────────────
SECRET_KEY = os.getenv('JWT_SECRET', 'fallback_secret_key') 

@app.route('/register', methods=['POST'])
def register():
    data = request.json
    name = data.get('name')
    email = data.get('email')
    password = data.get('password')
    
    if not name or not email or not password:
        return jsonify({'error': 'Missing fields'}), 400
        
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('SELECT id FROM users WHERE email = ?', (email,))
    if cursor.fetchone():
        conn.close()
        return jsonify({'error': 'You have already registered with this account!'}), 400
        
    password_hash = bcrypt.hashpw(password.encode('utf-8'), bcrypt.gensalt()).decode('utf-8')
    
    cursor.execute('INSERT INTO users (name, email, password_hash) VALUES (?, ?, ?)', (name, email, password_hash))
    conn.commit()
    conn.close()
    
    return jsonify({'message': 'Registration successful'}), 200

@app.route('/login', methods=['POST'])
def login():
    data = request.json
    email = data.get('email')
    password = data.get('password')
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('SELECT id, name, password_hash FROM users WHERE email = ?', (email,))
    user = cursor.fetchone()
    conn.close()
    
    if user and bcrypt.checkpw(password.encode('utf-8'), user[2].encode('utf-8')):
        token = jwt.encode({
            'user_id': user[0],
            'name': user[1],
            'exp': datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=24)
        }, SECRET_KEY, algorithm='HS256')
        return jsonify({'token': token, 'name': user[1]}), 200
        
    return jsonify({'error': 'Invalid credentials.'}), 401

@app.route('/forgot-password', methods=['POST'])
def forgot_password():
    data = request.json
    email = data.get('email')
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('SELECT id FROM users WHERE email = ?', (email,))
    if not cursor.fetchone():
        conn.close()
        return jsonify({'error': 'Email not registered.'}), 404
        
    otp = str(random.randint(100000, 999999))
    cursor.execute('UPDATE users SET reset_token = ? WHERE email = ?', (otp, email))
    conn.commit()
    conn.close()
    
    # Normally emit an email. We just return it for testing purposes.
    return jsonify({'message': 'OTP generated.', 'otp': otp}), 200

@app.route('/reset-password', methods=['POST'])
def reset_password():
    data = request.json
    email = data.get('email')
    otp = data.get('token')
    new_password = data.get('new_password')
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('SELECT id FROM users WHERE email = ? AND reset_token = ?', (email, otp))
    user = cursor.fetchone()
    
    if not user:
        conn.close()
        return jsonify({'error': 'Invalid OTP or email.'}), 400
        
    password_hash = bcrypt.hashpw(new_password.encode('utf-8'), bcrypt.gensalt()).decode('utf-8')
    cursor.execute('UPDATE users SET password_hash = ?, reset_token = NULL WHERE email = ?', (password_hash, email))
    conn.commit()
    conn.close()
    
    return jsonify({'message': 'Password reset successful.'}), 200

# ── Routes ───────────────────────────────────────────────────────────────────

@app.route('/')
def home():
    """Main landing page."""
    display_genes = [g for g in feature_cols if not g.endswith('_hits')]
    gene_info = {g: GENE_SIGS.get(g, {}) for g in display_genes}
    return render_template('index.html',
                           genes=display_genes,
                           gene_info=gene_info,
                           antibiotics=antibiotic_cols,
                           importances=IMPORTANCES)


@app.route('/predict', methods=['POST'])
def predict():
    """Predicts antibiotic resistance based on detected or selected genes."""
    try:
        if request.is_json:
            data = request.get_json()
        else:
            data = request.form

        screened_genes = BASE_FEATURE_COLS
        feature_subset = data.get('feature_subset', '')
        if feature_subset:
            try:
                requested = json.loads(feature_subset) if isinstance(feature_subset, str) else feature_subset
            except json.JSONDecodeError:
                return jsonify({'error': 'feature_subset must be valid JSON'}), 400

            if not isinstance(requested, list):
                return jsonify({'error': 'feature_subset must be a JSON array of gene names'}), 400

            screened_genes = [gene for gene in requested if gene in BASE_FEATURE_COLS]
            if not screened_genes:
                return jsonify({'error': 'feature_subset did not include any known genes'}), 400

        # 1. Get feature values (0 or 1)
        detected = {}
        for gene in feature_cols:
            if gene.endswith('_hits'):
                # Heuristically infer hits from binary gene presence if not provided
                base_gene = gene.rsplit('_hits', 1)[0]
                detected[gene] = int(data.get(gene, data.get(base_gene, 0)))
            else:
                detected[gene] = int(data.get(gene, 0))

        prediction_rows = predict_resistance(bundle, detected)

        results = []
        db_results = {}

        for row in prediction_rows:
            results.append({
                'antibiotic': row['antibiotic'],
                'prediction': row['prediction'],
                'confidence': row['confidence']
            })
            db_results[row['antibiotic']] = row['prediction']

        save_history = str(data.get('save_history', '1')).lower() in {'1', 'true', 'yes', 'on'}
        if save_history:
            filename = data.get('filename', 'Manual Input')
            conn = sqlite3.connect(DB_PATH)
            cursor = conn.cursor()
            col_names = [ab.replace('-', '_').replace(' ', '_') for ab in antibiotic_cols]
            placeholders = ", ".join(["?"] * (len(col_names) + 1))
            col_list = "filename, " + ", ".join(col_names)
            values = [filename] + [db_results.get(ab, 'N/A') for ab in antibiotic_cols]
            cursor.execute(
                f'INSERT INTO predictions ({col_list}) VALUES ({placeholders})',
                values
            )
            conn.commit()
            conn.close()

        return jsonify({
            'results': results,
            'genes_detected': [gene for gene, value in detected.items() if value == 1],
            'screened_genes': screened_genes,
            'saved_to_history': save_history
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/analyze-genome', methods=['POST'])
def analyze_genome():
    """Automatic gene detection from genome text/file."""
    try:
        raw = ''
        if 'fasta_file' in request.files:
            f = request.files['fasta_file']
            if f and f.filename:
                raw = f.read().decode('utf-8', errors='ignore')
        if not raw:
            raw = request.form.get('fasta', '').strip()
        
        if not raw:
            return jsonify({'error': 'No sequence provided.'}), 400

        feature_subset = request.form.get('feature_subset', '').strip()
        genes_to_screen = BASE_FEATURE_COLS
        if feature_subset:
            try:
                requested = json.loads(feature_subset)
            except json.JSONDecodeError:
                return jsonify({'error': 'Invalid feature subset payload.'}), 400

            if not isinstance(requested, list):
                return jsonify({'error': 'Feature subset must be a list of genes.'}), 400

            genes_to_screen = [gene for gene in requested if gene in BASE_FEATURE_COLS]
            if not genes_to_screen:
                return jsonify({'error': 'No valid genes were requested for screening.'}), 400

        dna = parse_fasta(raw)
        if len(dna) < 50:
            return jsonify({'error': 'Sequence too short for heuristic screening.'}), 400

        detected, meta = screen_resistance_markers(dna, genes_to_screen, GENE_SIGS)

        return jsonify({
            'detected': detected,
            'meta': meta,
            'genes_found': sum(1 for gene in genes_to_screen if detected.get(gene) == 1),
            'sequence_length': len(dna),
            'method': 'heuristic_reference_screen',
            'screened_genes': genes_to_screen,
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/history')
def history():
    """Returns past results."""
    try:
        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM predictions ORDER BY timestamp DESC')
        rows = cursor.fetchall()
        conn.close()
        return jsonify([dict(row) for row in rows])
    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/example-fasta')
def example_fasta():
    """Generates a sample genome with specific genes for demonstration."""
    requested = request.args.get('genes', 'blaNDM,gyrA_mut,tetA').split(',')
    lines = ['>E_coli_AMR_example | MicrobeNet Predictor | Synthetic genome']
    filler = 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC'
    seq = filler
    for gene in feature_cols:
        sigs = GENE_SIGS.get(gene, {}).get('signatures', [])
        seq += (sigs[0] if (gene in requested and sigs) else '') + filler
    for i in range(0, len(seq), 60):
        lines.append(seq[i:i+60])
    return jsonify({'fasta': '\n'.join(lines)})


# ── Smart Input Endpoints ────────────────────────────────────────────────────

@app.route('/api/get-disease-suggestions', methods=['GET'])
def get_disease_suggestions():
    """Returns a list of all available diseases for autocomplete."""
    query = request.args.get('q', '').lower()
    diseases = list(SMART_DATA.get('diseases', {}).keys())
    
    if query:
        matches = [d for d in diseases if query in d.lower()]
    else:
        matches = diseases
        
    return jsonify(matches)

@app.route('/api/get-bacteria-from-disease', methods=['GET'])
def get_bacteria_from_disease():
    """Returns a list of bacteria associated with a specific disease."""
    disease = request.args.get('disease', '')
    bacteria_list = SMART_DATA.get('diseases', {}).get(disease, [])
    return jsonify(bacteria_list)

@app.route('/api/get-sequence', methods=['GET'])
def get_sequence():
    """Returns the simulated genomic sequence for a specific bacteria."""
    bacteria_name = request.args.get('bacteria', '')

    matched_name, sequence, candidates = _find_bacteria_match(bacteria_name)

    if matched_name and sequence:
        return jsonify({
            'bacteria': matched_name,
            'sequence': _enrich_smart_sequence(matched_name, sequence),
            'markers': _resolve_profile_for_bacteria(matched_name),
            'screened_genes': _resolve_profile_for_bacteria(matched_name),
        })

    if candidates:
        return jsonify({
            'error': 'Ambiguous bacteria name. Please select a specific species.',
            'candidates': candidates,
        }), 409

    return jsonify({'error': 'Bacteria sequence not found'}), 404


@app.route('/api/get-bacteria-suggestions', methods=['GET'])
def get_bacteria_suggestions():
    """Returns a list of all available bacteria for autocomplete."""
    query = request.args.get('q', '').lower()
    bacteria_list = list(SMART_DATA.get('bacteria', {}).keys())
    
    if query:
        matches = [b for b in bacteria_list if query in b.lower()]
    else:
        matches = bacteria_list
        
    return jsonify(matches)



@app.route('/api/consult-docverse', methods=['POST'])
def consult_docverse():
    """Consult AI Docverse for alternative medicine suggestions."""
    data = request.json
    antibiotic = data.get('antibiotic')
    bacteria = data.get('bacteria', 'the target pathogen')
    
    if not antibiotic:
        return jsonify({'error': 'Antibiotic name is required'}), 400
        
    prompt = f"""As a clinical AI assistant named Docverse, provide alternative antibiotic suggestions or alternative medical treatments for a patient where the bacteria ({bacteria}) is resistant to {antibiotic}. 
    
    Structure your response with:
    1. A brief explanation of the resistance.
    2. 3-4 Alternative Antibiotics that are likely to be effective.
    3. Clinical considerations for these alternatives.
    
    Keep the tone professional, concise, and informative. Use Markdown formatting."""

    try:
        headers = {
            "Authorization": f"Bearer {GROQ_API_KEY}",
            "Content-Type": "application/json"
        }
        payload = {
            "model": GROQ_MODEL,
            "messages": [
                {"role": "system", "content": "You are Docverse, a professional clinical AI assistant specializing in antimicrobial resistance and alternative treatments. Provide expert medical insights based on antibiotic resistance profiles."},
                {"role": "user", "content": prompt}
            ],
            "temperature": 0.5,
            "max_tokens": 1024
        }
        
        response = requests.post(GROQ_BASE_URL, headers=headers, json=payload, timeout=20)
        response.raise_for_status()
        
        ai_response = response.json()
        content = ai_response['choices'][0]['message']['content']
        
        return jsonify({'suggestion': content})
        
    except requests.exceptions.HTTPError as http_err:
        print(f"Groq API Error Response: {http_err.response.text}")
        return jsonify({'error': f'AI Docverse API Error: {http_err.response.text}'}), 400
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Docverse AI Error Details: {str(e)}")
        return jsonify({'error': f'AI Docverse Error: {str(e)}'}), 500


if __name__ == '__main__':
    app.run(debug=False, port=5000)
