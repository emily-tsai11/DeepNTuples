import os

# ---------------------------------------------------
# CHANGE THESE CONFIG OPTIONS
submit_jobs = True
samples = ["TTToHadronic_PU200", "TTToHadronic_noPU"]
# ---------------------------------------------------

for sample in samples:
    print("Creating submit file for " + sample)

    submit_text  = "universe    = vanilla\n"
    submit_text += "executable  = run_job.sh\n"
    submit_text += "arguments   = $(inputFile) " + sample + " $(Process)\n"
    submit_text += "output      = status/" + sample + "_$(Process).out\n"
    submit_text += "log         = status/" + sample + "_$(Process).log\n"
    submit_text += "error       = status/" + sample + "_$(Process).err\n"
    submit_text += "+JobFlavour = \"longlunch\"\n"
    submit_text += "\n"
    submit_text += "queue inputFile from " + sample + ".list\n"

    submit_file = open(sample + ".submit", "w")
    submit_file.write(submit_text)
    submit_file.close()

    if submit_jobs:
        os.system("condor_submit " + sample + ".submit -batch-name " + sample)

if not submit_jobs:
    print("Jobs unsubmitted")
