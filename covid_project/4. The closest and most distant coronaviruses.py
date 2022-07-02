# Посчитать матрицу качеств попарных глобальных выравниваний геномов семи
# коронавирусов человека: SARS-CoV-2, SARS-CoV, MERS-CoV, HCoV-OC43, HCoV-NL63,
# HCoV-HKU1, HCoV-229E (fasta файлы можно найти в RefSeq).
# Найти пару ближайших и пару самых далеких вирусов.
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import itertools
import math

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 5
aligner.mismatch_score = -4
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

cov_lst = ['SARS-CoV-2', 'SARS-CoV', 'MERS-CoV', 'HCoV-OC43',
           'HCoV-NL63', 'HCoV-HKU1', 'HCoV-229E']
cov_dict = dict()
for elem in cov_lst:
    link = elem + '.fasta'
    cov_dict[elem] = list(SeqIO.parse(link, 'fasta'))[0].seq

scores = dict()
max_score = -math.inf
min_score = math.inf
max_score_pair = None
min_score_pair = None

all_pairs = list(itertools.combinations(cov_dict.keys(), 2))
for element in all_pairs:
    alignments = aligner.align(cov_dict[element[0]], cov_dict[element[1]])
    sc = float(alignments.score)
    scores[element] = sc
    if sc > max_score:
        max_score = sc
        max_score_pair = element
    if sc < min_score:
        min_score = sc
        min_score_pair = element

for e in scores:
    print(f'{e}: {scores[e]}')

print(f'Пара самых близких вирусов: {max_score_pair}, score: {max_score}\n'
      f'Пара самых далеких вирусов: {min_score_pair}, score: {min_score}')

# Вывод программы (21 парное выравнивание обрабатывается достаточно долго, a
# Jupiter выдает Memory Error, хотя в PyCharm все работает):
# ('SARS-CoV-2', 'SARS-CoV'): 95872.0
# ('SARS-CoV-2', 'MERS-CoV'): 41925.5
# ('SARS-CoV-2', 'HCoV-OC43'): 39885.0
# ('SARS-CoV-2', 'HCoV-NL63'): 36771.0
# ('SARS-CoV-2', 'HCoV-HKU1'): 40439.0
# ('SARS-CoV-2', 'HCoV-229E'): 35839.5
# ('SARS-CoV', 'MERS-CoV'): 41186.0
# ('SARS-CoV', 'HCoV-OC43'): 38636.5
# ('SARS-CoV', 'HCoV-NL63'): 35382.0
# ('SARS-CoV', 'HCoV-HKU1'): 37983.0
# ('SARS-CoV', 'HCoV-229E'): 34811.5
# ('MERS-CoV', 'HCoV-OC43'): 39944.0
# ('MERS-CoV', 'HCoV-NL63'): 35933.5
# ('MERS-CoV', 'HCoV-HKU1'): 40009.5
# ('MERS-CoV', 'HCoV-229E'): 34817.0
# ('HCoV-OC43', 'HCoV-NL63'): 38424.5
# ('HCoV-OC43', 'HCoV-HKU1'): 77791.0
# ('HCoV-OC43', 'HCoV-229E'): 36234.5
# ('HCoV-NL63', 'HCoV-HKU1'): 40779.5
# ('HCoV-NL63', 'HCoV-229E'): 60647.5
# ('HCoV-HKU1', 'HCoV-229E'): 36934.0
# Пара самых близких вирусов: ('SARS-CoV-2', 'SARS-CoV'), score: 95872.0
# Пара самых далеких вирусов: ('SARS-CoV', 'HCoV-229E'), score: 34811.5
