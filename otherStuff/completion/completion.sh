_testc_completions()
{
  COMPREPLY+=($(compgen -W  "$(ls /scratch/shared/egige60/ElmerIce_NEMO_RealCPL/PARAMETERS/)"  "${COMP_WORDS[1]}"))
}

complete -F _testc_completions ./createRun.sh

