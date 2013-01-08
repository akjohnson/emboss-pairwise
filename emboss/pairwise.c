/* @source pairwise application
**
** Pairwise comparisons of a set of sequences with
** True Needleman-Wunsch global alignment
**
** @author Audra Johnson (audrakjohnson@gmail.com)
** based on "needle" by
** @author Copyright (C) Alan Bleasby (ableasby@hgmp.mrc.ac.uk)
** @@
**
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
******************************************************************************/

#include "emboss.h"

/* @prog needle ***************************************************************
**
** Needleman-Wunsch global alignment
**
******************************************************************************/

int main(int argc, char **argv)
{
    AjPAlign align;
    AjPSeqset seqall;
    const AjPSeq a;
    const AjPSeq b;
    AjPStr m;
    AjPStr n;
    AjPStr ss;

    /* These help define sections of the pairwise alignment to do,
       allowing the work to be split up over a cluster more easily.
    */
    ajint start;
    ajint end;

    ajuint lena;
    ajuint lenb;

    const char *p;
    const char *q;

    ajint start1 = 0;
    ajint start2 = 0;

    float *path;
    ajint *compass;

    AjPMatrixf matrix;
    AjPSeqCvt cvt = 0;
    float **sub;

    float gapopen;
    float gapextend;
    ajulong maxarr = 1000;     /* arbitrary. realloc'd if needed */
    ajulong len;

    float score;

    AjBool dobrief = ajTrue;

    float id   = 0.;
    float sim  = 0.;
    float idx  = 0.;
    float simx = 0.;

    AjPStr tmpstr = NULL;

    size_t stlen;

    embInit("pairwise", argc, argv);

    matrix    = ajAcdGetMatrixf("datafile");
    seqall    = ajAcdGetSeqset("asequence");
    start     = ajAcdGetInt("start");
    end       = ajAcdGetInt("end");
    gapopen   = ajAcdGetFloat("gapopen");
    gapextend = ajAcdGetFloat("gapextend");
    dobrief   = ajAcdGetBool("brief");

    ajSeqsetTrim(seqall);

    align     = ajAcdGetAlign("outfile");

    gapopen = ajRoundF(gapopen, 8);
    gapextend = ajRoundF(gapextend, 8);

    AJCNEW(path, maxarr);
    AJCNEW(compass, maxarr);

    m  = ajStrNew();
    n  = ajStrNew();
    ss = ajStrNew();

    sub = ajMatrixfArray(matrix);
    cvt = ajMatrixfCvt(matrix);

    ajuint i;
    ajuint j;
    ajuint nseqs = ajSeqsetGetSize(seqall);

    //If no end specified, go to the end of the file
    if( end == -1 ) {
        end = nseqs;
    }
    else if( end > nseqs ) {
        ajDie( "The sequence at %d specified to end the pairwise comparisons, but only %d sequences in the file.\n", end, nseqs );
    }

    for( i = start - 1; i < end; i++ ) {
        a = ajSeqsetGetseqSeq(seqall, i);
        lena = ajSeqGetLen(a);

        for( j = i + 1; j < nseqs; j++ ) {
            b = ajSeqsetGetseqSeq(seqall, j);
            lenb = ajSeqGetLen(b);

            if( lenb > ( ULONG_MAX / (ajulong)(lena+1) ) )
                 ajFatal("Sequences too big. Try 'stretcher' or 'supermatcher'");

            len = lena * lenb;

            if( len > maxarr ) {
                stlen = (size_t) len;
                AJCRESIZETRY(path, stlen);

                if( !path )
                    ajDie("Sequences too big. Try 'stretcher'");

                AJCRESIZETRY(compass, stlen);

                if( !compass )
                    ajDie("Sequences too big. Try 'stretcher'");

                maxarr = len;
            }

            p = ajSeqGetSeqC(a);
            q = ajSeqGetSeqC(b);

            ajStrAssignC(&m, "");
            ajStrAssignC(&n, "");

            embAlignPathCalc(p, q, lena, lenb, gapopen, gapextend, path, sub, cvt,
                    compass, ajFalse);

            score = embAlignScoreNWMatrix(path, a, b, sub, cvt, lena, lenb,
                        gapopen, compass,
                        gapextend, &start1, &start2);

            embAlignWalkNWMatrix(path, a, b, &m, &n, lena, lenb, &start1, &start2, gapopen,
                            gapextend, cvt, compass, sub);

            embAlignReportGlobal(align, a, b ,m, n,
                             start1, start2,
                             gapopen, gapextend,
                             score, matrix,
                             ajSeqGetOffset(a), ajSeqGetOffset(b));

            if( !dobrief )
            {
                embAlignCalcSimilarity(m, n, sub, cvt, lena, lenb, &id, &sim, &idx,
                         &simx);
                ajFmtPrintAppS(&tmpstr,"Longest_Identity = %5.2f%%\n",
                     id);
                ajFmtPrintAppS(&tmpstr,"Longest_Similarity = %5.2f%%\n",
                     sim);
                ajFmtPrintAppS(&tmpstr,"Shortest_Identity = %5.2f%%\n",
                     idx);
                ajFmtPrintAppS(&tmpstr,"Shortest_Similarity = %5.2f%%",
                     simx);
                ajAlignSetSubHeaderApp(align, tmpstr);
            }
            ajAlignWrite(align);
            ajAlignReset(align);
        }
    }

    ajAlignClose(align);
    ajAlignDel(&align);

    ajSeqsetDel(&seqall);

    AJFREE(compass);
    AJFREE(path);

    ajStrDel(&n);
    ajStrDel(&m);
    ajStrDel(&ss);
    ajStrDel(&tmpstr);

    embExit();

    return 0;
}
