# This is a resubmission

## Comments

All the comments from the reviewer were addressed (email from Sept 5th):

> Thanks,
>
> Please do not start the description with "This package", package name,
> title or similar.

Fixed.

>
> If there are references describing the methods in your package, please
> add these in the description field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> authors (year, ISBN:...)
> or if those are not available: <https:...>
> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
> auto-linking.
> (If you want to add a title as well please put it in quotes: "Title")
>

Added some references (with DOIs) to the Description as suggested.


> Please always explain all acronyms in the description text. -> FFT

Fixed.

> Please always add all authors, contributors and copyright holders in the
> Authors@R field with the appropriate roles.
>  From CRAN policies you agreed to:
> "The ownership of copyright and intellectual property rights of all
> components of the package must be clear and unambiguous (including from
> the authors specification in the DESCRIPTION file). Where code is copied
> (or derived) from the work of others (including from R itself), care
> must be taken that any copyright/license statements are preserved and
> authorship is not misrepresented.
> Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the
> authors of such code. Alternatively, the ‘Author’ field should list
> these authors as contributors.
> Where copyrights are held by an entity other than the package authors,
> this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’
> field, or using a ‘Copyright’ field (if necessary referring to an
> inst/COPYRIGHTS file)."
> e.g.:
> Please explain in the submission comments what you did about this issue.
> 

Fixed. I explicitly mentioned all authors as contributors / copyright holders. This
package includes code derived from the R-base stats package. On my first
submission I listed "R Core Team and contributors worldwide" in `Authors@R`,
the CRAN reviewer asked me to list all authors and contributors explicitly.
Therefore I checked the copyright header of the files I used and the svn history
and I included all their full names with a role of both contributors and copyright
holders.

> Please fix and resubmit.
> 

Thanks for your detailed feedback.

## Test environments
* local Ubuntu 22.04 install, R 4.2.1
* win-builder (devel and release)
* r_hub

## R CMD check results
There are no NOTEs, WARNINGs or ERRORs.

## Downstream dependencies

This is a new package and there are no downstream dependencies.
