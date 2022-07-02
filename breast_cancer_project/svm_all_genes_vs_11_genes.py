import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, RocCurveDisplay

df = pd.read_pickle("bc_data.pkl")
an = pd.read_pickle("bc_ann.pkl")

with open('genes_SVM.txt') as f:
    genes = f.readlines()

genes = genes[0].strip('\n').split(';')
df_genes = df[genes]

X_train1 = df.loc[an["Dataset type"] == "Training"].to_numpy()
y_train1 = an.loc[an["Dataset type"] == "Training", "Class"].to_numpy()

X_test1 = df.loc[an["Dataset type"] == "Validation"].to_numpy()
y_test1 = an.loc[an["Dataset type"] == "Validation", "Class"].to_numpy()

X_train2 = df_genes.loc[an["Dataset type"] == "Training"].to_numpy()
y_train2 = an.loc[an["Dataset type"] == "Training", "Class"].to_numpy()

X_test2 = df_genes.loc[an["Dataset type"] == "Validation"].to_numpy()
y_test2 = an.loc[an["Dataset type"] == "Validation", "Class"].to_numpy()

model1 = SVC(kernel="linear", class_weight="balanced").fit(X_train1, y_train1)
y_pred1 = model1.predict(X_test1)

model2 = SVC(kernel="linear", class_weight="balanced").fit(X_train2, y_train2)
y_pred2 = model2.predict(X_test2)

model3 = SVC(kernel="linear", class_weight="balanced").fit(X_train1, y_train1)
y_pred3 = model1.predict(X_train1)

model4 = SVC(kernel="linear", class_weight="balanced").fit(X_train2, y_train2)
y_pred4 = model2.predict(X_train2)


def classification_quality(model, X_test, y_test, y_pred, name):
    M = confusion_matrix(y_test, y_pred)
    TP = M[1, 1]
    TN = M[0, 0]
    FN = M[1, 0]
    FP = M[0, 1]
    TPR = TP / (TP + FN)
    TNR = TN / (TN + FP)
    RocCurveDisplay.from_estimator(model, X_test, y_test)
    plt.savefig(f"{name}.png", dpi=300)
    print(f'TPR {name}: {TPR}\n'
          f'TNR {name}: {TNR}\n')


classification_quality(model1, X_test1, y_test1, y_pred1, 'all_genes_model_test')
classification_quality(model2, X_test2, y_test2, y_pred2, '11_genes_model_test')
classification_quality(model3, X_train1, y_train1, y_pred3, 'all_genes_model_train')
classification_quality(model4, X_train2, y_train2, y_pred4, '11_genes_model_train')

# TPR all_genes_model_test: 0.4583333333333333
# TNR all_genes_model_test: 0.75
#
# TPR 11_genes_model_test: 0.7083333333333334
# TNR 11_genes_model_test: 0.8055555555555556
#
# TPR all_genes_model_train: 1.0
# TNR all_genes_model_train: 1.0
#
# TPR 11_genes_model_train: 0.7301587301587301
# TNR 11_genes_model_train: 0.7912087912087912

# при построении модели на выборке генов, наблюдаем лучшее качество по всем
# показателям. возможно, выбранные гены в большей степени связаны с рецидивом рака
# груди, а остальные гены, наоборот, добавляют шум в данные, из-за чего страдает
# качество классификации при построении модели на данных всех генов
# при обучении на всех генов и проверки качества на той же обучающей выборке
# классификатор получился идеальным. а при обучении только на части генов он
# показывает неплохое качество, однако видимо признаков недостаточно для идеальной классификации
