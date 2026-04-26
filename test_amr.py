import os
import sys
sys.path.append('.')
from amr_utils import load_model_bundle, predict_resistance, screen_resistance_markers

def test_prediction():
    try:
        bundle = load_model_bundle('.')
        feature_cols = bundle['feature_cols']
        print(f"Total features: {len(feature_cols)}")
        
        # Test case 1: All zeros
        detected_0 = {g: 0 for g in feature_cols}
        results_0 = predict_resistance(bundle, detected_0)
        print("\nResults for all zeros:")
        for r in results_0[:3]:
            print(f"  {r['antibiotic']}: {r['prediction']} ({r['confidence']}%)")
            
        # Test case 2: blaNDM = 1
        detected_ndm = {g: 0 for g in feature_cols}
        detected_ndm['ndm_beta_lactamase'] = 1
        if 'ndm_beta_lactamase_hits' in feature_cols:
            detected_ndm['ndm_beta_lactamase_hits'] = 5
            
        results_ndm = predict_resistance(bundle, detected_ndm)
        print("\nResults for blaNDM = 1:")
        for r in results_ndm[:3]:
            print(f"  {r['antibiotic']}: {r['prediction']} ({r['confidence']}%)")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_prediction()
