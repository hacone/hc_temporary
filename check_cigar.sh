#/bin/bash

gawk 'BEGIN{ FS="\t"; print "RID","RPOS","CIGAR","NM","MD" } $6~/^[0-9]*[HS][0-9]*M[0-9]*[HS]$/{ print $3,$4,$6,$7,$8 }' .tmp.mapped > .tmp.SMS.cigars

gawk 'BEGIN{ FS="\t"; print "RID","RPOS","CIGAR","NM","MD" } $6~/^[0-9]*M[0-9]*[HS]$/{ print $3,$4,$6,$7,$8 }' .tmp.mapped > .tmp.MS.cigars

gawk 'BEGIN{ FS="\t"; print "RID","RPOS","CIGAR","NM","MD" } $6~/^[0-9]*[HS][0-9]*M$/{ print $3,$4,$6,$7,$8 }' .tmp.mapped > .tmp.SM.cigars

gawk 'BEGIN{ FS="\t"; print "RID","RPOS","CIGAR","NM","MD" } $6~/^([0-9]*[HS])*([0-9]*M[0-9]*D)+[0-9]*M([0-9]*[HS])*$/{ print $3,$4,$6,$7,$8 }' .tmp.mapped > .tmp.D.cigars

gawk 'BEGIN{ FS="\t"; print "RID","RPOS","CIGAR","NM","MD" } $6~/^[0-9]*M$/{ print $3,$4,$6,$7,$8 }' .tmp.mapped > .tmp.M.cigars
