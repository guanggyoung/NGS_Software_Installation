
/******************************************************************
 * From a sorted TSV file containing bases found from sc-BS-Seq experiment,
 * merge the same base into one single line and write to anther two files:
 *  _merged.tsv and _merged.bed
 * 
 * Call the compiled program in this way
 * merged_sorted_tsv_keep_strand_info sample_concatenate_sorted.tsv
 * 
 * with get 2 files: 
 * sample_concatenate_sorted.tsv_merged.tsv
 * sample_concatenate_sorted.tsv_merged.bed
 * 
 * The program handle the TSV file line by line. And write the processed line instantly to the disk.
 * Because very big TSV (>4G) will use > 64GB RAM if read in the whole TSV directly int RAM.
 * 
 *****************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> //for stringstream
#include <algorithm> // for std::min() function

struct TSV_entry
{
	std::string ref_chr_name; // chromosome name
	char strand; // '+' or '-'
	unsigned long int position; // 1-based left-most mapping position
	
	float methylated_bases; 
	float un_methylated_bases;

	float met_ratio;

	// overloading "<<" operator to print SNP lines from TSV_file ifstream
	friend std::ostream & operator << (std::ostream & TSV_output,  TSV_entry & one_TSV_line)
	{
		TSV_output<<one_TSV_line.ref_chr_name<<'\t';
		TSV_output<<one_TSV_line.position<<'\t';
        TSV_output<<one_TSV_line.strand<<'\t';
		TSV_output<<one_TSV_line.methylated_bases<<'\t';
        TSV_output<<one_TSV_line.un_methylated_bases<<'\t';
        TSV_output<<one_TSV_line.met_ratio<<std::endl;

		return TSV_output;
	}
	
};

void merge_same_base_in_TSV(std::string TSV_filename)
{
       
    std::ifstream TSV_input(TSV_filename);
	std::string one_string_line;
	
	// std::vector<TSV_entry> whole_TSV; // NOT used
	
    std::getline(TSV_input, one_string_line); // get the first line but don't accept the info: Table head line
    
    // Write the new TSV file with sequence info into a new file.
    std::ofstream updated_TSV(TSV_filename+"_merged.tsv");
    // Write the new bed file with sequence info into a new file.
    std::ofstream updated_bed(TSV_filename+"_merged.bed");
    // print table head
    updated_TSV<<"ref"<<'\t'<<"pos"<<'\t'<<"strand"<<'\t'<<"methylated"<<'\t'<<"un_methylated"<<'\t'<<"met_ratio"<<std::endl;
    updated_bed<<"ref"<<'\t'<<"start"<<'\t'<<"end"<<'\t'<<"strand"<<'\t'<<"met_ratio"<<std::endl;
    
	std::getline(TSV_input, one_string_line);
	std::stringstream ss(one_string_line);
    TSV_entry old_TSV_line; //Use this to get the first 
	ss>>old_TSV_line.ref_chr_name;
	ss>>old_TSV_line.strand;
	ss>>old_TSV_line.position;
    ss>>old_TSV_line.methylated_bases;
    ss>>old_TSV_line.un_methylated_bases;
    ss>>old_TSV_line.met_ratio;



    // begin opration on the TSV file line by line
	while( std::getline(TSV_input, one_string_line) )
	{
		std::stringstream ss(one_string_line);
        TSV_entry one_TSV_line; // define one_TSV_line to get the info in this line
		ss>>one_TSV_line.ref_chr_name;
		ss>>one_TSV_line.strand;
		ss>>one_TSV_line.position;
        
        ss>>one_TSV_line.methylated_bases;
        ss>>one_TSV_line.un_methylated_bases;
        ss>>one_TSV_line.met_ratio;

        if(one_TSV_line.ref_chr_name == old_TSV_line.ref_chr_name &&
           one_TSV_line.position == old_TSV_line.position &&
           one_TSV_line.strand == old_TSV_line.strand) // same base in the genome, merge the new line to old
        {
            old_TSV_line.methylated_bases = old_TSV_line.methylated_bases + one_TSV_line.methylated_bases;
            old_TSV_line.un_methylated_bases = old_TSV_line.un_methylated_bases + one_TSV_line.un_methylated_bases;
        }
        else // A different base found, write the old base into updated 
        {   
            old_TSV_line.met_ratio = old_TSV_line.methylated_bases/(old_TSV_line.methylated_bases + old_TSV_line.un_methylated_bases);
            updated_TSV<<old_TSV_line;
            updated_bed<<old_TSV_line.ref_chr_name<<'\t'<<old_TSV_line.position<<'\t'<<old_TSV_line.position<<'\t'
                       <<old_TSV_line.strand<<'\t'<<old_TSV_line.met_ratio<<std::endl;

            old_TSV_line = one_TSV_line;
        }
    }
    updated_TSV.close();
    updated_bed.close();
}
/***************************************************************
argv[1] get TSV filename
call the program in this way:
.merged_sorted_tsv_keep_strand_info sample_concatenate_sorted.tsv
****************************************************************/

int main(int argc, char* argv[])
{
    std::string TSV_filename = argv[1];

    merge_same_base_in_TSV(TSV_filename);

    return 0;
}
