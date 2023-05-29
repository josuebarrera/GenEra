#!/bin/bash

cd test_files
mkdir tmp_test

gunzip test_DB.dmnd.gz

genEra -q test_seqs.fasta \
 -t 559292 -b test_DB \
 -d test_taxdump \
 -r test_ncbi_lineages.csv \
 -n 1 -x tmp_test \
 -s test_distances.tsv -i true > log

WARNINGTEST=0

GENERASCRIPT=$(which genEra)
if [[ -z ${GENERASCRIPT} ]]; then
	echo
	echo "  ERROR: Please make sure that genEra is located in your PATH"
	echo
	gzip test_DB.dmnd
	rm -fr tmp_test
	cd ../
	exit 1
fi

PHYLOSCRIPT=$(which Erassignment)
if [[ -z ${PHYLOSCRIPT} ]]; then
	echo
	echo "  ERROR: Please make sure that Erassignment is located in your PATH"
	echo
	gzip test_DB.dmnd
	rm -fr tmp_test
	cd ../
	exit 1
fi


if grep -q "YBR238C	XP_033764827.1	0.0	1350	27291" tmp_test/tmp_559292_*/559292_Diamond_results.bout; then
	echo "STEP 1: PASSED"
else
	echo "STEP 1: FAILED"
	echo "please make sure that diamond is working correctly"
	echo "or contact josue.barrera@tuebingen.mpg.de to address the issue"
	gzip test_DB.dmnd
	rm -fr tmp_test
	cd ../
	exit 1
fi

if grep -q "1071382,Kazachstania africana CBS 2517,Kazachstania,Saccharomycetaceae,Saccharomycetales,saccharomyceta" 559292_ncbi_lineages.csv; then
	echo "STEP 2: PASSED"
else
	echo "STEP 2: FAILED"
	echo "Please make sure to follow the github instructions for a correct installation"
	echo "or contact josue.barrera@tuebingen.mpg.de to address the issue"
	gzip test_DB.dmnd
	rm -fr tmp_test
	cd ../
	exit 1
fi

if grep -q "2	Saccharomycetales	2" 559292_gene_age_summary.tsv; then
	echo "STEP 3: PASSED"
else
	echo "STEP 3: FAILED"
	echo "Please make sure to follow the github instructions for a correct installation"
	echo "or contact josue.barrera@tuebingen.mpg.de to address the issue"
	gzip test_DB.dmnd
	rm -fr tmp_test
	cd ../
	exit 1
fi

if grep -q "YBR238C,YGL107C	Saccharomycetales	2	2" 559292_founder_events.tsv; then
	echo "gene clustering: PASSED"
else
	echo "gene clustering: FAILED"
	echo "please make sure that mcl is working correctly"
	echo "or contact josue.barrera@tuebingen.mpg.de to address the issue"
	gzip test_DB.dmnd
	rm -fr tmp_test
	cd ../
	exit 1
fi

if grep -q "YBR238C	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	Ortholog_detected	0.0" 559292_abSENSE_results/Detection_failure_probabilities; then
	echo "STEP 4: PASSED"
else
	echo "STEP 4: FAILED"
	echo "WARNING: Homology detection failure cannot be tested"
	echo "Please make sure that abSENSE.py is working correctly"
	echo "or contact josue.barrera@tuebingen.mpg.de to address the issue"
	WARNINGTEST+=1
fi

rm -fr tmp_test 559292_* log
gzip test_DB.dmnd

if [ ${WARNINGTEST} -eq 0 ]; then
	echo
	echo "GenEra is ready to be used!"
	echo
else
	echo
	echo "Found ${WARNINGTEST} WARNINGS, genEra results may be suboptimal"
	echo
fi

cd ../

exit 0
