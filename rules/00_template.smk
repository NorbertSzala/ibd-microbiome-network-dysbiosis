rule rulename:
    
    """
    Template rule

    Description:
        Srutututu
    
    """
    input:

    output:
    
    params:

    conda:
        "../envs/.yaml"

    log:
        f"{LOGS}/template.smk"

    message:
        "Runnign template rule"

    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname ${output.XXX}) $(dirname {log})


        run script \
            -I \
            -O \
            > {log} 2>&1


        """
