from IPython.utils import process
from numpy.random.mtrand import random_integers
import numpy as np
from sklearn.linear_model import LogisticRegression
import torch
import random
import math
import warnings
from sklearn import preprocessing
#from google.colab import drive
import sys
#drive.mount('/content/drive')
#import os
#os.chdir('drive/MyDrive')
warnings.filterwarnings(action='ignore')

class CustomDataset(torch.utils.data.Dataset):
    def __init__(self, statements, truth):
        # TODO: initialize self.statements and self.truth
        self.statements = statements
        self.truth = truth
        
    def __getitem__(self, index):
        # TODO: return self.statement and self.truth at the given index as a tuple
        return self.statements[index], self.truth[index]
    
    def __len__(self):
        # TODO: return the number of elements in dataset
        return len(self.truth)
bounded = open('bound.fasta')
unbounded = open('notbound.fasta')
tests = open("test.fasta")

def separate_columns(strings):
    columns = []

    # Find the maximum number of columns
    max_columns = max(len(string) for string in strings)

    for i in range(max_columns):
        column = []

        for string in strings:
            if i < len(string):
                column.append(string[i])
            else:
                column.append('')

        columns.append(column)

    return columns



def getText(file, txt_type):
  statements =[]
  labels =[]
  text = []
  index = 0
  countBound = 0
  countUnbound = 0
  if txt_type == "bound":
    txt = ">bound"
  elif txt_type == "unbound":
    txt = ">notbound"
  else:
    txt = ">seq"
  for i in file:
    if txt not in i:
      text.append(i[:len(i)-1])
    elif len(text) > 0:
      statements.append(text)
      text = []
      if txt_type == "bound":
        labels.append(1)
      else:
        labels.append(0)
      index+=1
  return statements, labels

statements, labels = getText(bounded, "bound")
unboundStatements, unboundLabels = getText(unbounded, "unbound")
seqText = getText(tests, "seq")
seqText = seqText[0]
statements.append(unboundStatements)
labels.append(unboundLabels)
for i in range(len(statements)): #Shuffle dataset
  rand = random.randint(0,len(statements)-1)
  temp = statements[i]
  tempLabel = labels[i]
  temp
  statements[i] = statements[rand]
  labels[i] = labels[rand]
  labels[rand] = tempLabel
  statements[rand] = temp

def create_pwm(sequences):
    # Count the occurrences of each nucleotide at each position
    pwm = {
        'A': [],
        'C': [],
        'G': [],
        'T': []
    }

    sequence_length = len(sequences[0])
    num_sequences = len(sequences)+4

    for i in range(sequence_length):
        position_counts = {
            'A': 1,
            'C': 1,
            'G': 1,
            'T': 1
        }

        for sequence in sequences:
            nucleotide = sequence[i]
            position_counts[nucleotide] += 1

        # Calculate the frequencies and add them to the PWM
        for nucleotide in ['A', 'C', 'G', 'T']:
            frequency = position_counts[nucleotide] / num_sequences
            pwm[nucleotide].append(frequency)

    return pwm



    '''

pwmInput = [{}]
for i in statements:
  for k in i:
    for j in range(len(statements)):
      pwmInput[j][k] += 1
      '''
eachCol = []
allCols = []
pwm =[]
for i in statements:
  sepCol = separate_columns(i)
  eachCol.append(sepCol)

for i in eachCol:
  if pwm != []:
    allCols.append(pwm)
    pwm = []
  for j in i:
    processedPWM=create_pwm(j)
    pwm.append(processedPWM)
cols = np.array(allCols)
labels = np.asarray(labels)

dataset = CustomDataset(statements,labels)
maxAcc = 0
bestRatio = 0
split_ratio = 0.3
train_dataset = dataset[:math.floor(len(dataset)*split_ratio)]
test_dataset = dataset[math.ceil(len(dataset)*split_ratio):]

train_loader = torch.utils.data.DataLoader(train_dataset, 128, shuffle=True)
test_loader = torch.utils.data.DataLoader(test_dataset, 128, shuffle=True)



'''
labelEnc = preprocessing.LabelEncoder()
new_target = labelEnc.fit_transform(statements)
onehotEnc = preprocessing.OneHotEncoder()
onehotEnc.fit(new_target.reshape(-1, 1))
targets_trans = onehotEnc.transform(new_target.reshape(-1, 1))
statements = targets_trans.toarray()

X_train = statements[:math.floor(len(statements)*split_ratio)].reshape(-1,1)
y_train = labels[:math.floor(len(labels)*split_ratio)]
X_test = statements[math.ceil(len(statements)*split_ratio):].reshape(-1,1)
y_test = labels[math.ceil(len(labels)*split_ratio):]
# Create and train the logistic regression model
model = LogisticRegression(max_iter=10000)
model.fit(X_train, y_train)

# Make predictions on the testing set
predictions = model.predict(X_test)
count  = 0
for i in y_test:
  if i == 0:
    count += 1
print(count/len(y_test))
for i in predictions:
  if i != 0:
    print(i)
# Evaluate the model
accuracy = np.mean(predictions == y_test)
print("Accuracy: "+str(accuracy))

bound_sequences = X_test[predictions == 1]  # Assuming '1' represents the bound class
 #Print the bound DNA sequences
for sequence in bound_sequences:
  finalSeq = ""
  for i in sequence:
    match i:
      case 1:
        finalSeq += 'A'
      case 2:
        finalSeq += 'C'
      case 3: 
        finalSeq += 'G'
      case 4:
        finalSeq += 'T'
  print(finalSeq)
'''
def create_profile_matrix(sequences):
    base_counts = {'A': [], 'C': [], 'G': [], 'T': []}
    profile_matrix = {}

    sequence_length = len(sequences[0])
    sequence_count = len(sequences) +4

    for base in base_counts:
        for i in range(sequence_length):
            base_counts[base].append(1)

    for sequence in sequences:
        for i in range(sequence_length):
            base = sequence[i]
            base_counts[base][i] += 1

    for base in base_counts:
        profile_matrix[base] = [count / sequence_count for count in base_counts[base]]

    return profile_matrix


def find_profile_most_probable_kmer(profile_matrix, sequence):
    k = 21
    max_probability = -1
    profile_most_probable_kmer = ""

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        probability = 1.0

        for j in range(k):
            base = kmer[j]
            probability *= profile_matrix[base][j]

        if probability > max_probability:
            max_probability = probability
            profile_most_probable_kmer = kmer

    return profile_most_probable_kmer

def compute_probability(profile_matrix, kmer):
    probability = 1.0

    for j in range(len(kmer)):
        base = kmer[j]
        probability *= profile_matrix[base][j]

    return probability

def rank_sequences(profile_matrix, sequences, top_n):
    probabilities = []

    for sequence in sequences:
        kmer = find_profile_most_probable_kmer(profile_matrix, sequence)
        probability = compute_probability(profile_matrix, kmer)
        probabilities.append((sequence, probability))

    probabilities.sort(key=lambda x: x[1], reverse=True)
    top_sequences = probabilities[:top_n]

    return top_sequences


floatTxt = ''
row = []
mat = []
import re
pwm = open("project2a_PWM.txt")
nFound = False
for i in pwm:
  for j in i:

    if str(j) != ' ' and str(j) != '\n' and str(j) != '\t':
      if str(j) == "." and "." in floatTxt :
            float1 = floatTxt[:len(floatTxt)-1]
            row.append(float(float1))
            float1 = ""
            floatTxt = "0."
            continue
      floatTxt += str(j)
    elif (j == '\t' or j == '\n') and floatTxt != "":
      row.append(float(floatTxt))
      floatTxt = ""
      if j == '\n':
        nFound = True
    elif  nFound and len(row) > 0:
      mat.append(row)
      row = []
      nFound = False
mat.append(row)

matDict = {}
arr = ["A", "C", "G", "T"]
for j in range(len(mat)):
    matDict[arr[j]] = mat[j]

'''

adenine = []
cytosine = []
guanine = []
thymine = []
for i in dictArr:
  for key in i:
    if key == 'A':
      adenine.append(i[key])
    if key == 'C':
      cytosine.append(i[key])
    if key == 'G':
      guanine.append(i[key])
    if key == 'T':
      thymine.append(i[key])
      finalDict = {'A': adenine, 'C': cytosine, 'G': guanine, 'T': thymine}

      '''

profile_matrix = []
for i in statements:
  cols = separate_columns(i)
  for j in cols:
    profile_matrix.append(create_profile_matrix(j))

profile_matrix = np.asarray(profile_matrix)

top_sequences = []
for i in range(len(seqText)):
    top_sequences.append(rank_sequences(matDict, seqText[i], 1))

seqDic = {}
for i in top_sequences:
  if "[" not in str(i[0][0]):
    seqDic[str(i[0][0])] = i[0][1]

sortedDic = dict(sorted(seqDic.items(), key=lambda item: item[1], reverse=True))
top_seqs = []
for key in sortedDic:
  top_seqs.append(key)
top_seqs = top_seqs[:2000]

for i in top_seqs:
  print(i)
