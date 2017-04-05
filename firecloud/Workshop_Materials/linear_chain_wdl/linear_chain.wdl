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

task multiplyByFiveTask {

	String inputNum

	command <<<
		OUT_NUM=$((${inputNum} * 5)) ;
		echo $OUT_NUM
		>>>

	output {
		String outNum=read_string(stdout())
		}
	
	runtime {
		docker: "ubuntu:14.04.4"
		}

	}


task addThreeTask {

	String inputNum

	command <<<
		OUT_NUM=$((${inputNum} + 3)) ;
		echo $OUT_NUM 
		>>>

	output {
		String final=read_string(stdout())
		}
	
	runtime {
		docker: "ubuntu:14.04.4"
		}

	}

workflow calculatorWorkflow {

	String inputNum

	call addTwoTask {
		input:
			inputNum=inputNum
		}

	call multiplyByFiveTask {
		input:
			inputNum=addTwoTask.outNum
		}

	call addThreeTask {
		input:
			inputNum=multiplyByFiveTask.outNum
		}


	}
