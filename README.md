# ige-ege

Project code and data from the article [A Bayesian decision-making model of implicit motor learning from internal and external errors](https://www.biorxiv.org/content/10.1101/2025.01.30.635749v1).


1) To install the same Python packages and versions I used for my analyses, you can recreate the same virtual environment by entering the following command in your terminal (within your cloned/downloaded project directory):

    `$ conda env create --name <nameyourenvironment> --file environment.yml`

2) These notebooks call on a package of custom-written functions found in the `src` directory. You may need to pip install this package. You can do so with the following command:

    `pip install -e .`


3) Once you've installed and conda activated the replicated environment, and pip installed the source code, you should be able to seamlessly run the Jupyter notebooks found in `experiment-1/scripts` to recreate my analyses:
    - `analyze-behavior.ipynb`
    - `model-fitting.ipynb`
    - `parameter-model-recovery.ipynb`
      

4) The necessary data and results for running some parts of these notebooks are found in `experiment-1/results/`.


Correspondence: hyosub.kim [at] ubc [dot] ca









