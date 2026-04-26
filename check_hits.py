import joblib
import os
try:
    cols = joblib.load('model/feature_cols.pkl')
    hits_cols = [c for c in cols if c.endswith('_hits')]
    print("Total features:", len(cols))
    print("Hits features count:", len(hits_cols))
    print("First 10 hits features:", hits_cols[:10])
except Exception as e:
    print("Error:", e)
