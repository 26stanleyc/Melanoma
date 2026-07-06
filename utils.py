import pickle

with open('/Users/stanleychen/git/Melanoma/MultiCoFusion/data/processed_data.pkl', 'rb') as f:
    data_cv = pickle.load(f)

print("Top-level keys:", data_cv.keys())
print("cv_splits keys:", data_cv['cv_splits'].keys())
print("Train Samples:", data_cv['cv_splits']['train'][:5])