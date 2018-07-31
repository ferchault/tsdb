import sys
import random
import cPickle
import numpy as np
import common

if __name__ == "__main__":
  parameters = {
  "train" : "train.txt",
  "test" : "test.txt",
	"kfold" : True;
  "N" : [10,20,40,80],
  "CV" : 10,
  "lambda" : 1e-7,
  "sigmas" : [0.1*2**i for i in range(25)],
  "Representation" : "SLATM",
  "Kernel" : "Laplacian",
  "seeed" : 667
  }

  random.seed(parameters["seeed"])
  nSigmas = len(parameters["sigmas"])

  # get properties
  data = common.get_properties(parameters["property_file"])
  mols, total = common.get_molecules(sorted(data.keys()), data)

  # get representaion
  X, Y = common.get_Representation(mols, parameters["Representation"])
  print parameters["sigmas"][2]

  # FCHL Kernel (sigmas implemented)
  '''if parameters["Representation"] == "FCHL":
    K = common.get_Kernel(X, parameters["sigmas"], parameters["Kernel"])
    for j in range(nSigmas):
      maes = []
      for train in parameters["N"]:
        maes, s = comcom.CrossValidation_fchl(K, train, parameters["CV"], parameters["lambda"], j, total, Y)
        s = np.std(maes)/np.sqrt(parameters["CV"])
        print(str(parameters["sigmas"][j]) +  "\t" + str(train) + "\t" + str(sum(maes)/len(maes)) + " " + str(s))
  # all other Kernels (sigma seperate)
  else:
  '''

  for j in range(nSigmas):
    K = common.get_Kernel(X, parameters["sigmas"], parameters["Kernel"])
    maes = []

    for train in parameters["N"]:
      maes, s = common.CrossValidation(K, train, parameters["CV"], parameters["lambda"], total, Y, parameters["kfold"])
      s = np.std(maes)/np.sqrt(parameters["CV"])

      print(str(parameters["sigmas"][j]) + "\t" + str(train) + "\t" + str(sum(maes)/len(maes)) + " " + str(s))

  # save Kernel
  print "\n save kernel"
  with open("K.cpickle", 'wb') as f:
    cPickle.dump(K, f, protocol=2)
