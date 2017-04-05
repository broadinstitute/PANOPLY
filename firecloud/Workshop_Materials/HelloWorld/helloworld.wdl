
workflow helloWorldWorkflow {
	String name_WF
	call helloWorldTask
	{
		input: name_T=name_WF
	}
}

task helloWorldTask {
	String name_T

command <<<
# Note: 
#    (1) must provide full pathname because working directory is not the root of the file system
#    (2) python and hello_world.py are bundled into the docker image cbirger/hello-world:5.0
python /hello_world/hello_world.py ${name_T}
>>>

output {
	File Hello_World_OutputFile = "hello_world_output.txt"
}

runtime {
	docker: "cbirger/hello-world:5.0"
}

}
