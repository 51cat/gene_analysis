import os
import argparse
import subprocess

from add_log import add_log

from collections import Counter

from cobra.io import read_sbml_model
from cobra.io import load_json_model
from cobra.io import load_matlab_model
from cobra.flux_analysis import single_gene_deletion

@add_log
def load_model_from_file(file : str, **kwargs):
    '''
    Parameters
    file: Your model file path

    Returns
    A cobra.model objective

    Raise
    FileNotFoundError
    IOError
    '''
    if not os.path.exists(file):
        raise FileNotFoundError(f"Not Found {file}")
    
    file_ext = file.lower().split('.')[-1]
    
    load_model_from_file.logger.info(f"{file_ext}")
    
    funcs_dict = {
        "xml": read_sbml_model,
        "json": load_json_model,
        "mat": load_matlab_model
    }
    
    try:
        model_obj = funcs_dict[file_ext](file, **kwargs)
    except KeyError:
        raise FileNotFoundError(f"error {file_ext}")
    except IOError:
        raise IOError(f"error {file}")
    
    return model_obj


class Knockout:
    def __init__(
            self, 
            model_file: str, 
            outdir:     str,
            method:     str = "fba") -> None:
        '''
        ## Parameters
        model_file:     Your model file path
        outdir:         Output directory
        method:         Method used to predict the growth rate {"fba", "moma", "linear moma", "room", "linear room"}
        
        ## Cmdline Usage
        python path/to/Knockout/knockout.py
                --model_file path/to/model_file
                --outdir path/to/output_directory
                --method fba default = fba
                --gene_list path/to/gene_list.txt

        ## Main Output
        {model_file_name}_knock_res.txt

        ## A tsv File
        ids	growth	status	name
        PP_1614	0.0	optimal	ispD
        PP_1167	0.5	optimal	PP_1167
        PP_3163	0.8	optimal	benC

        '''    
        self.model_file = model_file
        self.outdir = outdir
        self.method = method

        self.gene_list = None
        self.model_obj = None
        self.gene_info = {}
        self.out_name = self.model_file.split("/")[-1].split(".")[0]
        self.out_res = f"{self.outdir}/{self.out_name}_knock_res.txt"

        if not os.path.exists(self.outdir):
            subprocess.check_call(f"mkdir -p {self.outdir}", shell=True)

    
    def load_data(self):
        self.model_obj = load_model_from_file(self.model_file)
        self.parse_gene_info_from_model()
    
    def set_gene_list(self, gene_list: list):
        '''
        gene_list: A list of genes 
        '''
        self.gene_list = gene_list
    
    @add_log
    def parse_gene_info_from_model(self):
        for r in self.model_obj.genes:
            self.gene_info.update({r.id:r.name})
        self.parse_gene_info_from_model.logger.info(f"Total Genes: {len(self.gene_info)}")

    @add_log
    def run_knockout(self):
        gene_status = single_gene_deletion(
                model = self.model_obj,
                gene_list = self.gene_list,
                method = self.method
        )
        gene_status["ids"] = gene_status["ids"].apply(lambda x: list(x)[0])
        gene_status["name"] = gene_status["ids"].apply(lambda x: self.gene_info[x])
        self.run_knockout.logger.info(f"{Counter(gene_status['status'].to_list())}")
        gene_status.to_csv(self.out_res, sep = "\t", index=None)


def main():
    parser = argparse.ArgumentParser(description='Single gene deletion analysis')
    parser.add_argument('--model_file', type=str, help='Path of model file', required=True)
    parser.add_argument('--outdir', type=str, help='Output directory', required=True)
    parser.add_argument('--method', type=str, 
                        help='Method used to predict the growth rate {"fba", "moma", "linear moma", "room", "linear room"}', default="fba")
    parser.add_argument('--gene_list', type=str, help='Each line is a gene_id must be within the model', default=None)
    args = parser.parse_args() 
    
    runner = Knockout(
        model_file=args.model_file,
        outdir=args.outdir,
        method=args.method
    )

    if args.gene_list is not None:
        with open(args.gene_list) as fd:
            gene_list = [gene.strip("\n") for gene in fd.readlines()]
            runner.set_gene_list(gene_list)
    
    runner.load_data()
    runner.run_knockout()


if __name__ == '__main__':
    main()