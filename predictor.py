#!/usr/bin/env python3
"""
MicrobeNet AMR Predictor — CLI Tool
Usage:
    python predictor.py <genome.fasta>
    python predictor.py <genome.fasta> --ollama
    python predictor.py <genome.fasta> --genes blaNDM mecA gyrA_mut
"""
import argparse
import json
import os
import sys
import urllib.error
import urllib.request

from amr_utils import load_model_bundle, parse_fasta, predict_resistance, screen_resistance_markers

BASE = os.path.dirname(os.path.abspath(__file__))

STATUS = {'Resistant': '❌ Resistant', 'Intermediate': '⚠️  Intermediate', 'Susceptible': '✅ Susceptible'}
COLORS = {'Resistant': '\033[91m', 'Intermediate': '\033[93m', 'Susceptible': '\033[92m'}
RESET = '\033[0m'
BOLD = '\033[1m'
DIM = '\033[2m'


def _extract_ollama_json_payload(text: str) -> dict:
    start = text.find('{')
    end = text.rfind('}') + 1
    if start == -1 or end <= start:
        raise RuntimeError('Ollama response did not contain a JSON object')

    try:
        return json.loads(text[start:end])
    except json.JSONDecodeError as exc:
        raise RuntimeError(f'Ollama returned malformed JSON: {exc}') from exc


def detect_by_ollama(dna: str, feature_cols: list, gene_sigs: dict,
                     ollama_url: str = 'http://localhost:11434') -> tuple[dict, str]:
    req = urllib.request.Request(f'{ollama_url}/api/tags')
    resp = urllib.request.urlopen(req, timeout=3)
    data = json.loads(resp.read())
    models = data.get('models', [])
    if not models:
        raise RuntimeError('No Ollama models installed. Run: ollama pull llama3.2')
    model_name = models[0]['name']

    gene_desc = '\n'.join(
        f'- {g}: {gene_sigs.get(g, {}).get("family", g)}'
        for g in feature_cols
    )
    gene_list = ', '.join(feature_cols)
    json_template = '{' + ', '.join(f'"{g}": 0' for g in feature_cols) + '}'
    prompt = (
        f"You are an expert clinical microbiologist. "
        f"Identify which AMR genes are present in this genome.\n\n"
        f"Genes: {gene_list}\n\nDescriptions:\n{gene_desc}\n\n"
        f"FASTA (first 3000 bp):\n{dna[:3000]}\n\n"
        f"Reply ONLY with a JSON object: {json_template}"
        f"\nUse 1=present, 0=absent. No explanations."
    )
    body = json.dumps({'model': model_name, 'prompt': prompt,
                       'stream': False, 'options': {'temperature': 0.1}}).encode()
    req = urllib.request.Request(
        f'{ollama_url}/api/generate', data=body,
        headers={'Content-Type': 'application/json'}, method='POST'
    )
    resp = urllib.request.urlopen(req, timeout=90)
    result = json.loads(resp.read())
    text = result.get('response', '')
    parsed = _extract_ollama_json_payload(text)
    return {g: int(bool(parsed.get(g, 0))) for g in feature_cols}, model_name


def print_results(results: list, genes_detected: list, importances: dict, method: str):
    ab_width = max(len(r['antibiotic']) for r in results) + 2

    print(f'\n{BOLD}  MicrobeNet AMR Predictor{RESET}')
    print(f'  {"─" * 56}')

    present = [g for g in genes_detected]
    if present:
        print(f'  {DIM}Genes detected:{RESET} {", ".join(present)}')
    else:
        print(f'  {DIM}Genes detected:{RESET} None')
    print(f'  {DIM}Method: {method}{RESET}')
    print(f'\n  {"─" * 56}')
    print(f'  {BOLD}{"Antibiotic":<{ab_width}} {"Prediction":<16} {"Confidence":>10}{RESET}')
    print(f'  {"─" * 56}')

    for row in results:
        pred = row['prediction']
        conf = row['confidence']
        color = COLORS.get(pred, '')
        label = STATUS.get(pred, pred)
        bar = '█' * int(conf / 10) + '░' * (10 - int(conf / 10))
        print(f'  {row["antibiotic"]:<{ab_width}} {color}{label:<16}{RESET} {conf:>5.1f}%  {DIM}{bar}{RESET}')

    if importances:
        top = sorted(importances.items(), key=lambda item: -item[1])[:5]
        print(f'\n  {"─" * 56}')
        print(f'  {BOLD}Top Resistance Drivers{RESET}')
        for gene, importance in top:
            bar = '█' * int(importance * 50)
            pct = importance * 100
            print(f'  {gene:<14} {pct:5.1f}%  {DIM}{bar}{RESET}')

    print(f'  {"─" * 56}\n')


def main():
    parser = argparse.ArgumentParser(
        description='MicrobeNet AMR Predictor — Predict antibiotic resistance from FASTA genome'
    )
    parser.add_argument('fasta', help='Path to FASTA genome file')
    parser.add_argument('--ollama', action='store_true',
                        help='Use Ollama AI for gene detection (requires Ollama running)')
    parser.add_argument('--genes', nargs='+', metavar='GENE',
                        help='Manually specify present genes (skips detection)')
    args = parser.parse_args()

    try:
        bundle = load_model_bundle(BASE)
    except FileNotFoundError:
        print('❌  Model not found. Run: python scripts/train_per_antibiotic_models.py')
        sys.exit(1)
    feature_cols = bundle['feature_cols']
    gene_sigs = bundle['gene_defs']
    importances = bundle['importances']

    try:
        with open(args.fasta, encoding='utf-8') as handle:
            raw = handle.read()
    except FileNotFoundError:
        print(f'❌  File not found: {args.fasta}')
        sys.exit(1)

    dna = parse_fasta(raw)
    print(f'\n  📂 File: {args.fasta}  ({len(dna):,} bp)')

    method = 'manual'
    if args.genes:
        detected = {}
        for g in feature_cols:
            if g.endswith('_hits'):
                base_gene = g.rsplit('_hits', 1)[0]
                detected[g] = int(base_gene in args.genes)
            else:
                detected[g] = int(g in args.genes)
    elif args.ollama:
        print('  🤖 Querying Ollama AI...')
        try:
            detected, model_name = detect_by_ollama(dna, feature_cols, gene_sigs)
            method = f'ollama:{model_name}'
        except Exception as exc:
            print(f'  ⚠️  Ollama failed ({exc}) → falling back to heuristic screening')
            detected, _ = screen_resistance_markers(dna, feature_cols, gene_sigs)
            method = 'heuristic_fallback'
    else:
        detected, _ = screen_resistance_markers(dna, feature_cols, gene_sigs)
        method = 'heuristic_reference_screen'

    results = predict_resistance(bundle, detected)
    genes_detected = [g for g, v in detected.items() if v == 1]
    print_results(results, genes_detected, importances, method)


if __name__ == '__main__':
    main()
