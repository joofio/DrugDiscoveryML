# DrugDiscoveryML
messing around with drug discovery and machine-learning

### Purpose
The objective of this repository is to address pharmacokinetics issues with drug discovery and if/how can machine-learning help on this issue while learning myself. Any input is welcomed.

### How to
The import cleans the data and imports a subset from chembl_23 and eache one of the assays works from there.
The dataset used is available [here](https://www.dropbox.com/s/jmhxpdn9m3izt02/product_adme.csv?dl=0).

### Assays
There are testing issues for logP, Protein Binding and Aqueous solubility

### Results so far
#### Solubility 

The mean is -2.145955  
The median is -2.25  

| Model     | RMSE     | MAE     | info                      | 
|-----------|----------|---------|---------------------------|
| glm stack | 0.959    | 0.524   | kknn + cubist + bam + gam |


#### LogP 

The mean of the y is 5.503537  
The median is 5  
![alt text](/images/logp.png "LogP results")

#### Binding energy

The mean is 10.11311  
The median is 11.2  
![alt text](/images/Binding.png "Binding Energy results")

### TODO
* Bioavailability
* Clean code and readability
* Improve results
  * Scale and normalize
  * Better feature selection
  * More methods
* ~~Introduce DrugBank Information~~
* Finishing test in python with scikit-learn
* Introduce KNIME, deepchem and RDkit and CNN
* ...
