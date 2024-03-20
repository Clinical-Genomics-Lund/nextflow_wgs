<!--
Thanks for contributing to the CMD nextflow_wgs pipeline!

Please use the checklists below to document changes and performed tests.

Remember that this template doubles as our test/review documentation. 
-->
# Description

What is changed? What is the result of the changes? 
How does this update improve the pipeline?

## Risk assessment

What are the possible challenges and issues that might result from this update? How likely and severe are they? 

## Type of change
<!--
    Major change counts as a change that breaks backward compatibility
    Minor change is a substantial change that requires testing before deployment
    Patch is a minor change like a bug fix, code comment/style fix, etc.
    
    Choose one and delete the remaining fields.
-->
- [ ] Documentation
- [ ] Patch
- [ ] Minor change
- [ ] Major change 

# Checklist:
<!--
    The checklist below applies to all types of changes, 
    except documentation updates. 

    In the case of changes that only affect documentation external to 
    code source files, remove the checkboxes related to code review
    and testing for no new warnings.
        
    Do not hesitate to add your own items to the checklist if applicable.
-->
- [ ] My code follows the [style guidelines for NF pipelines at CMD](http://mtlucmds1.lund.skane.se/wiki/doku.php?id=nextflow&s[]=nextflow#code_style_at_cmd)
- [ ] I have performed a self-review of my code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] My changes generate no new warnings (Keep an eye on `.nextflow.log` !)
- [ ] I have updated the CHANGELOG
- [ ] The latest commit in the master branch is tagged
      with the correct version number

<!--
    Select a checklist below based on selection under # Type of change
    and delete the sections that do not apply to this PR:
-->
## Major change
- [ ] Stub run completes without errors or new warnings
- [ ] `onco` run finishes without any new warnings/errors and the results can 
       be loaded into scout
- [ ] `wgs` single run finishes without any new warnings/errors and the results 
       can be loaded into scout
- [ ] `wgs` trio run finishes without any new warnings/errors and the results 
       can be loaded into scout
- [ ] At least one other person has reviewed and approved my code
- [ ] I have made corresponding changes to the documentation

## Minor change
- [ ] Stub run completes without errors or new warnings
- [ ] `onco` run finishes without any new warnings/errors and the results can 
       be loaded into scout
- [ ] `wgs` single run finishes without any new warnings/errors and the results 
      can be loaded into scout
- [ ] `wgs` trio run finishes without any new warnings/errors and the results 
       can be loaded into scout
- [ ] At least one other person has reviewed and approved my code
- [ ] I have made corresponding changes to the documentation

## Patch
- [ ] Stub run completes without errors or new warnings
<!--
- [ ] At least one other person has reviewed and approved my code
-->

<!--
## [Optional] Are there any post-deployment tasks we need to perform?

 - [ ] Task/link to issue 1
 - [ ] Task/link to issue 2
 - [ ] ...
-->

# Instructions for the reviewers

## How to test the changes
<!--
    Provide clear, concise and specific steps for reviewers to test your changes.    
-->
1. Step 1
2. Step 2

# Test/review documentation

## Review performed by:

- [ ] Alexander
- [ ] Jakob
- [ ] Paul
- [ ] Ryan
- [ ] Viktor
    
## Testing performed by:

- [ ] Alexander
- [ ] Jakob
- [ ] Paul
- [ ] Ryan
- [ ] Viktor
