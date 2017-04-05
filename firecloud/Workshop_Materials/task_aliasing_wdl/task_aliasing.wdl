task addTwoTask {

	String inputNum

	command <<<
		OUT_NUM=$((${inputNum} + 2)) ;
		echo $OUT_NUM
		>>>

	output {
		String outNum=read_string(stdout())
		}

	runtime {
		docker: "ubuntu:14.04.4"
		}
	}

task addNumsTask {

	String inputNumOne
	String inputNumTwo

	command <<<
		OUT_NUM=$((${inputNumOne} + ${inputNumTwo})) ;
		echo $OUT_NUM > outFile.txt
		echo "The numbers $IN_NUM_ONE and $IN_NUM_TWO were added!" > log.txt
		>>>

	output {
		File outFile="outFile.txt"
		String outNum=read_string("outFile.txt") 
		File log="log.txt"
		}
	
	runtime {
		docker: "ubuntu:14.04.4"
		}
	}

workflow TaskAliasingWorkflow {

	String dataOne
	String dataTwo

	call addTwoTask as FirstAdder {
		input:
			inputNum=dataOne
		}


	call addTwoTask as SecondAdder {
		input:
			inputNum=dataTwo
		}

	call addNumsTask {
		input:
			inputNumOne=FirstAdder.outNum,
			inputNumTwo=SecondAdder.outNum
		}

}


