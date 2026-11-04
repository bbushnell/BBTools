package prok;

import java.util.Arrays;
import java.util.Random;

import map.IntHashMap;
import structures.IntList;

/** TaxID-grouped validation split protecting rare phyla and original-model bins.
 * Input contains EVERY original locus, before partial/excluded learning filters.
 * Each TaxID is indivisible and belongs to one phylum. Validation targets one
 * tenth of groups, subject to the same cap in each phylum and model bin.
 * Arrays supplied by the caller are borrowed and never modified.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningSplit{

	/** Returns a nonempty feasible validation partition; underfilling the target is allowed. */
	public static Result split(final int[] taxid, final int[] phylum, final int[] model, final int phylumCount, final int modelCount, final long seed){
		if(taxid==null || phylum==null || model==null || taxid.length==0
				|| phylum.length!=taxid.length || model.length!=taxid.length || phylumCount<1 || modelCount<1){
			throw new IllegalArgumentException("Split requires aligned original-locus arrays and nonempty phylum/model axes");
		}
		final IntHashMap lookup=new IntHashMap();
		final IntList unique=new IntList();
		for(int gene=0; gene<taxid.length; gene++){
			if(taxid[gene]<=0 || phylum[gene]<0 || phylum[gene]>=phylumCount || model[gene]<0 || model[gene]>=modelCount){
				throw new IllegalArgumentException("Invalid original TaxID/phylum/model at locus "+gene);
			}
			final int previous=lookup.get(taxid[gene]);
			if(previous<0){
				lookup.put(taxid[gene], phylum[gene]);
				unique.add(taxid[gene]);
			}else if(previous!=phylum[gene]){
				throw new IllegalArgumentException("One TaxID cannot cross phyla: "+taxid[gene]);
			}
		}
		final int[] groups=unique.toArray();
		Arrays.sort(groups);
		final int target=groups.length/10;
		if(target==0){throw new IllegalArgumentException("Validation needs at least ten positive TaxID groups");}
		final int[] groupPhylum=new int[groups.length];
		final int[] phylumTotal=new int[phylumCount];
		for(int group=0; group<groups.length; group++){
			groupPhylum[group]=lookup.get(groups[group]);
			phylumTotal[groupPhylum[group]]++;
			// After collecting phylum assignments, this map becomes TaxID -> group index.
			lookup.put(groups[group], group);
		}
		final int words=(int)(((long)modelCount+63)/64);
		final long cells=(long)words*groups.length;
		if(cells>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Group/model membership exceeds indexed-array capacity");}
		final long[] membership=new long[(int)cells];
		final int[] modelTotal=new int[modelCount];
		final int[] groupByGene=new int[taxid.length];
		for(int gene=0; gene<taxid.length; gene++){
			final int group=lookup.get(taxid[gene]);
			groupByGene[gene]=group;
			final int cell=group*words+(model[gene]>>>6);
			final long bit=1L<<(model[gene]&63);
			if((membership[cell]&bit)==0){
				membership[cell]|=bit;
				modelTotal[model[gene]]++;
			}
		}
		final int[] order=new int[groups.length];
		for(int group=0; group<order.length; group++){order[group]=group;}
		final Random random=new Random(seed);
		for(int i=order.length-1; i>0; i--){
			final int j=random.nextInt(i+1);
			final int value=order[i];
			order[i]=order[j];
			order[j]=value;
		}
		final boolean[] validation=new boolean[groups.length];
		final int[] phylumValidation=new int[phylumCount];
		final int[] modelValidation=new int[modelCount];
		int selected=0;
		for(int group:order){
			if(selected==target){break;}
			final int p=groupPhylum[group];
			if(phylumValidation[p]>=phylumTotal[p]/10
					|| blocked(group, words, membership, modelTotal, modelValidation)){
				continue;
			}
			validation[group]=true;
			phylumValidation[p]++;
			selected++;
			for(int word=0; word<words; word++){
				long bits=membership[group*words+word];
				while(bits!=0){
					final int m=(word<<6)+Long.numberOfTrailingZeros(bits);
					modelValidation[m]++;
					bits&=bits-1;
				}
			}
		}
		if(selected==0){throw new IllegalArgumentException("No validation group fits the rare-phylum/model protection limits");}
		assert(selected<=target) : "Validation may underfill, never exceed, its group target";
		for(int p=0; p<phylumCount; p++){
			assert(phylumValidation[p]<=phylumTotal[p]/10) : "Validation removed more than one tenth of phylum "+p;
		}
		for(int m=0; m<modelCount; m++){
			assert(modelValidation[m]<=modelTotal[m]/10) : "Validation removed more than one tenth of original model bin "+m;
		}
		final boolean[] validationByGene=new boolean[taxid.length];
		for(int gene=0; gene<taxid.length; gene++){validationByGene[gene]=validation[groupByGene[gene]];}
		return new Result(groups, groupPhylum, groupByGene, validation, validationByGene, target, selected);
	}

	/** Checks only bins represented by this group; repeated copies never multiply group counts. */
	private static boolean blocked(final int group, final int words, final long[] membership, final int[] total, final int[] validation){
		assert(total.length==validation.length) : "Model quotas and validation counts must share the original model axis";
		for(int word=0; word<words; word++){
			long bits=membership[group*words+word];
			while(bits!=0){
				final int model=(word<<6)+Long.numberOfTrailingZeros(bits);
				if(validation[model]>=total[model]/10){return true;}
				bits&=bits-1;
			}
		}
		return false;
	}

	/** Owns its arrays; consumers treat them as read-only. TaxID groups are sorted. */
	public static final class Result{
		private Result(final int[] taxids_, final int[] phyla_, final int[] groupByGene_, final boolean[] validation_,
				final boolean[] validationByGene_, final int target_, final int selected_){
			taxids=taxids_;
			phylumByGroup=phyla_;
			groupByGene=groupByGene_;
			validationByGroup=validation_;
			validationByGene=validationByGene_;
			target=target_;
			selected=selected_;
		}
		public final int[] taxids, phylumByGroup, groupByGene;
		public final boolean[] validationByGroup, validationByGene;
		public final int target, selected;
	}

	private NcrnaLearningSplit(){}
}
