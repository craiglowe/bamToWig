/*

bamToWig.c

*/

#include "common.h"
#include "options.h"
#include "memalloc.h"
#include "hash.h"
#include "linefile.h"
#include "bed.h"
#include "sam.h"
#include "bamFile.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"bedGraph", OPTION_BOOLEAN},
	{"expansion", OPTION_INT},
	{NULL, 0}
};

int optExpansion = 0;
boolean optBedGraph = FALSE;

/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	"bamToWig - create a wiggle representing the read depth of the bam\n"
	"usage:\n"
	"   bamToWig noGap.bed in.bam output.wig\n"
	"options:\n"
	"   -expansion         (0)        Number of bases to extend from start site of read\n"
	"                                  instead of using the end of the read\n"
	);
}

/*---------------------------------------------------------------------------*/

struct slChrom
{
	struct slChrom *next;
	char *name;
	unsigned int length;
};


struct slChrom *createChrom(char *name, unsigned int length)
{
	struct slChrom *ret;
	AllocVar(ret);
	ret->name = cloneString(name);
	ret->length = length;
	return(ret);
}

struct hash *bedLoadNInHash(char *filename, int fields, char *restrictToChrom)
{
	struct bed *regions = NULL, *bed = NULL, *currList = NULL, *temp = NULL, *nextBed = NULL;
	struct hash *regionHash = newHash(6);

	regions = bedLoadNAll(filename, fields);
	slSort(&regions, bedCmp);
	currList = regions;
	for(bed = regions; bed != NULL; bed = nextBed)
	{
		nextBed = bed->next;
		if((bed->next == NULL) || (differentString(bed->chrom,bed->next->chrom)))
		{
			temp = bed->next;
			bed->next = NULL;
			if(restrictToChrom == NULL || sameString(bed->chrom, restrictToChrom))
			{
				hashAdd(regionHash, bed->chrom, currList);
			}
			else
			{
				bedFreeList(&currList);
			}
			currList = temp;
		}
	}
	return(regionHash);
}


struct slChrom *chromListFromNoGapHash(struct hash *noGapHash)
{
	struct hashCookie cookie = hashFirst(noGapHash);
	unsigned int max = 0;
	struct bed *headBed = NULL, *bunk = NULL;
	struct slChrom *headChrom = NULL, *currChrom = NULL;
	struct hashEl *hel = NULL;

	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		headBed = hel->val;
		max = 0;
		for(bunk = headBed; bunk != NULL; bunk=bunk->next)
		{
			if(bunk->chromEnd > max){max = bunk->chromEnd;}
		}
		currChrom = createChrom(headBed->chrom, max);
		currChrom->next = headChrom;
		headChrom = currChrom;
	}
	return(headChrom);
}


struct hash *chromListToUnsignedHash(struct slChrom *chromList)
{
	struct slChrom *currChrom = NULL;
	struct hash *answer = newHash(6);

	for(currChrom=chromList; currChrom != NULL; currChrom=currChrom->next)
	{
		unsigned int *baseData = NULL;
		AllocArray(baseData, currChrom->length);
		hashAddUnique(answer, currChrom->name, baseData);
	}

	return(answer);
}


unsigned int findChromMax(struct hash *noGapHash, char *chrom)
{
	struct bed *curr = NULL;
	unsigned int max = 0;
	for(curr=hashFindVal(noGapHash, chrom); curr != NULL; curr=curr->next)
	{
		if(curr->chromEnd > max){max=curr->chromEnd;}
	}
	return(max);
}


unsigned int minUnsigned(unsigned int a, unsigned int b)
{
	if(a <= b){return(a);}
	else{return(b);}
}


void addReadCounts(struct hash *noGapHash, struct hash *coverageHash, char *filename)
{
	samfile_t *samFp = NULL;
	samFp = samopen(filename, "rb", 0);
	bam1_t* b = bam_init1();
	bam1_core_t *c = NULL;
	char *chrom = NULL;
	int chromId = -1;
	unsigned int x=0, chromStart=0, chromEnd=0, chromMax=0;
	unsigned int *coverage = NULL;

	while(samread(samFp, b) >= 0)
	{
		c = &b->core;
		if(!(c->flag&(BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FUNMAP)))
		{
			if(c->tid != chromId)
			{
				chromId = c->tid;
				chrom = samFp->header->target_name[chromId];
				coverage = hashFindVal(coverageHash, chrom);
				chromMax = findChromMax(noGapHash, chrom);
			}

			chromStart = c->pos;
			chromEnd = bam_calend(c, bam1_cigar(b));
			if(optExpansion != 0)
			{
				if(c->flag&BAM_FREVERSE)
				{
					if(chromEnd > optExpansion){chromStart = chromEnd - optExpansion;}
					else{chromStart = 0;}
				}
				else
				{
					chromEnd = minUnsigned(chromStart + optExpansion, chromMax);
				}
			}
			if(coverage != NULL)
			{
				for(x = chromStart; x < chromEnd; x++)
				{
					coverage[x]++;
				}
			}
		}
	}
	samclose(samFp);
	bam_destroy1(b);
}


void printCoverage(struct hash *noGapHash, struct hash *coverageHash, char *outFilename)
{
	unsigned int basePos = 0;
	unsigned int *coverage = NULL;

	FILE *fout = mustOpen(outFilename, "w");

	struct hashCookie cookie = hashFirst(noGapHash);
	struct bed *headBed = NULL, *bunk = NULL;
	struct hashEl *hel = NULL;

	verbose(3, " Writing output...\n");
	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		headBed = hel->val;
		coverage = hashMustFindVal(coverageHash, headBed->chrom);
		for(bunk = headBed; bunk != NULL; bunk=bunk->next)
		{
			fprintf(fout, "fixedStep chrom=%s start=%u step=1\n", bunk->chrom, bunk->chromStart + 1);
			for(basePos=bunk->chromStart; basePos < bunk->chromEnd; basePos++)
			{
				fprintf(fout, "%u\n", coverage[basePos]);
			}
		}
	}
	carefulClose(&fout);
}


void printCoverageBedGraph(struct hash *noGapHash, struct hash *coverageHash, char *outFilename)
{
	unsigned int basePos = 0, printStart = 0, printNum = 0;
	unsigned int *coverage = NULL;

	FILE *fout = mustOpen(outFilename, "w");

	struct hashCookie cookie = hashFirst(noGapHash);
	struct bed *headBed = NULL, *bunk = NULL;
	struct hashEl *hel = NULL;

	verbose(3, " Writing output...\n");
	for(hel = hashNext(&cookie); hel != NULL; hel = hashNext(&cookie))
	{
		headBed = hel->val;
		coverage = hashMustFindVal(coverageHash, headBed->chrom);
		for(bunk = headBed; bunk != NULL; bunk=bunk->next)
		{
			printStart = bunk->chromStart;
			printNum = coverage[printStart];
			for(basePos=bunk->chromStart+1; basePos < bunk->chromEnd; basePos++)
			{
				if(coverage[basePos] != printNum)
				{
				fprintf(fout, "%s\t%u\t%u\t%u\n", bunk->chrom, printStart, basePos, printNum);
					printStart = basePos;
					printNum = coverage[basePos];
				}
				if(basePos == bunk->chromEnd-1)
				{
					fprintf(fout, "%s\t%u\t%u\t%u\n", bunk->chrom, printStart, basePos+1, printNum);
				}
			}
		}
	}
	carefulClose(&fout);
}


/*---------------------------------------------------------------------------*/


void bamToWig(char *noGapFilename, char *bamFilename, char *outFilename)
{
	struct hash *coverageHash = NULL, *noGapHash = NULL;
	struct slChrom *chromList = NULL;

	verbose(2, "Reading in no gap file\n");
	noGapHash = bedLoadNInHash(noGapFilename, 3, NULL);

	verbose(2, "Creating a list of chroms from the no gap file\n");
	chromList = chromListFromNoGapHash(noGapHash);

	verbose(2, "Creating data structure for mismap data\n");
	coverageHash = chromListToUnsignedHash(chromList);

	verbose(2, "Reading bam files\n");
	addReadCounts(noGapHash, coverageHash, bamFilename);

	verbose(2, "Writing output file\n");
	if(optBedGraph){printCoverageBedGraph(noGapHash, coverageHash, outFilename);}
	else{printCoverage(noGapHash, coverageHash, outFilename);}
}


/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 4){usage();}

	optExpansion = optionInt("expansion", optExpansion);
	optBedGraph = optionExists("bedGraph");

	bamToWig(argv[1], argv[2], argv[3]);

	return(0);
}

