import joblib
import os
try:
    cols = joblib.load('model/feature_cols.pkl')
    print("Features:", cols)
    print("Count:", len(cols))
except Exception as e:
    print("Error:", e)
