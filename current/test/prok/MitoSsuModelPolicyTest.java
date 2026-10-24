package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import map.LongHashSet;
import prok.*;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Exercises actual loading, short-window scheduling, verification and raw endpoints.
 * @author Brian, Raiden, Ganyu */
public final class MitoSsuModelPolicyTest {
	public static void main(final String[] args) throws Exception{
		check(args.length==1, "Supply a new fixture directory");final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		if(ml.Function.TYPE_RATES_CUM==null){ml.Function.normalizeTypeRates();}
		check(MitoSsuModelPolicy.minimumLength(50,.66f)==33 && MitoSsuModelPolicy.minimumLength(800,.66f)==528,
			"Float boundary support must not reject33/50 at the caller's own .66f threshold");
		final NcrnaFamily packaged=NcrnaFamilyLoader.load("16S");Prok16sReleaseTest.checkPackaged(packaged);
		// Project exactly the registry's original prefix, so synthetic additions never
		// inherit any production coverage model or silently weaken baseline checks.
		final NcrnaFamily original=new NcrnaFamily("16S",Arrays.copyOf(packaged.library,5),null,
			Arrays.copyOf(packaged.modelNames,5),Prok16sReleaseTest.originalSeeds(),17,1400,1800,9,1,false,0,0,0,1,0,1,.66f,.66f);
		final byte[][] library=Arrays.copyOf(original.library,10);final String[] names=Arrays.copyOf(original.modelNames,10);
		final Random random=new Random(72021);final int[] lengths={800,820,1935,790,810};
		for(int i=0; i<5; i++){library[i+5]=randomBases(random,lengths[i]);names[i+5]=MitoSsuModelPolicy.ADDITIONS[i];}
		final String path=dir.resolve("expanded.fa").toString();final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{for(int i=0; i<library.length; i++){out.print(b.clear().append('>').append(names[i]).nl().append(library[i]).nl());}}
		finally{RrnaResourceIO.finish(out);}
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();options.resources.consensus=path;
		final NcrnaFamily expanded=NcrnaFamilyLoader.load("16S",options,new FamilyResourceResolver());
		check(MitoSsuModelPolicy.expanded(expanded.name,expanded.modelNames) && expanded.library.length==10,
			"The actual unified loader must admit exactly the complete original plus five-model roster");
		for(int i=0; i<10; i++){
			final boolean raw=i==4 || i>=5;
			check(expanded.rrnaEndpointFeatures.raw[i]==raw && (expanded.rrnaEndpointFeatures.fiveNets[i]==null)==raw
				&& (expanded.rrnaEndpointFeatures.threeNets[i]==null)==raw,
				"The original four trained endpoint pairs stay active; only named mitochondrial models retain raw ends");
		}
		final LongHashSet seeds=new LongHashSet(8);seeds.add(key(library[5],80));seeds.add(key(library[5],350));
		final NcrnaScavenger caller=caller(library,names,seeds);caller.setRrnaEndpointFeatures(expanded.rrnaEndpointFeatures,null);
		for(int flank:new int[]{0,400}){
			final byte[] contig=new byte[800+2*flank];Arrays.fill(contig,(byte)'N');System.arraycopy(library[5],0,contig,flank,800);
			final ArrayList<Orf> calls=caller.scavenge("short_mito",contig,0,new ArrayList<int[]>());
			check(calls.size()==1 && calls.get(0).start==flank && calls.get(0).stop==flank+799
				&& calls.get(0).trnaModel.equals(names[5]) && calls.get(0).toGff().contains("model:mitochondria_chordate"),
				"A complete800-base mitochondrial gene must survive both whole-contig and interior windows, preserve raw ends and label the model");
		}
		check(caller.prepareCandidateWindows(new int[]{100,350},null,MitoSsuModelPolicy.minimumLength(790,.66f)).size()==1,
			"Collapsed windows at the inclusive mitochondrial floor must survive before model alignment");
		check(caller.scavenge("too_short",Arrays.copyOf(library[5],500),0,new ArrayList<int[]>()).isEmpty(),
			"The model-relative scheduling floor must still reject undersized contigs");
		final byte[][] oldLibrary=Arrays.copyOf(library,5);oldLibrary[0]=library[5];
		final NcrnaScavenger unchanged=caller(oldLibrary,Arrays.copyOf(names,5),seeds);
		check(unchanged.scavenge("old_floor",library[5],0,new ArrayList<int[]>()).isEmpty(),
			"An original-roster model must retain the existing prokaryotic input floor");
		final byte[] padded=new byte[1600];Arrays.fill(padded,(byte)'N');System.arraycopy(library[5],0,padded,400,800);
		check(unchanged.scavenge("old_alignment_floor",padded,0,new ArrayList<int[]>()).isEmpty(),
			"An800-base exact match to an old model must still fail its alignment-span floor inside a long window");
		final byte[][] mixed=library.clone();mixed[0]=randomBases(random,800);
		final LongHashSet oldSeeds=new LongHashSet(8);oldSeeds.add(key(mixed[0],80));oldSeeds.add(key(mixed[0],350));
		System.arraycopy(mixed[0],0,padded,400,800);
		check(caller(mixed,names,oldSeeds).scavenge("old_model_expanded_roster",padded,0,new ArrayList<int[]>()).isEmpty(),
			"With all ten models active, a perfect undersized original-model winner must still fail its original span floor");
		for(int length:new int[]{540,527}){
			final byte[] part=Arrays.copyOfRange(library[5],800-length,800);
			// Interior substitutions force the full-model identity below .66 while
			// keeping a high-identity physical suffix for the clipped verifier.
			for(int i=0; i<30; i++){final int pos=60+i*7;part[pos]=part[pos]=='A' ? (byte)'C' : (byte)'A';}
			final LongHashSet partSeeds=new LongHashSet(8);partSeeds.add(key(part,300));partSeeds.add(key(part,400));
			for(int strand=0; strand<2; strand++){
				final ArrayList<Orf> clipped=caller(library,names,partSeeds).scavenge("mito_partial",part,strand,new ArrayList<int[]>());
				if(length==527){check(clipped.isEmpty(),"A527-base clipped suffix must fail the800-base model's528-base support floor");}
				else{
					check(clipped.size()==1 && clipped.get(0).ncrnaPartial==1 && clipped.get(0).start==0
						&& clipped.get(0).stop==length-1 && clipped.get(0).trnaModel.equals(names[5]),
						"A540-base high-identity physical suffix must use genuine clipping with raw bounds and missing5prime provenance");
					final Orf call=clipped.get(0);if(strand==1){call.flip();}
					check(call.toGff().contains(";partial="+(strand==0 ? "10" : "01")),
						"Mitochondrial partial provenance must preserve genomic orientation in actual GFF formatting");
				}
			}
		}
		final String[] bad=names.clone();bad[9]="mitochondria_typo";
		check(!MitoSsuModelPolicy.expanded("16S",bad) && !MitoSsuModelPolicy.expanded("23S",names),
			"The policy cannot activate from an arbitrary prefix or in another family");
		final NcrnaScavenger incompatible=caller(library,names,seeds);incompatible.voteEnds=true;
		boolean rejected=false;try{incompatible.scavenge("invalid",library[5],0,new ArrayList<int[]>());}catch(IllegalStateException expected){rejected=true;}
		check(rejected,"Unsupported verification or voting must fail instead of silently sharing mitochondrial length policy");
		coverageTiers(original,dir,random);
		tierLimits(original,dir,random);
		System.out.println("MITO_SSU_POLICY_TEST_PASS loader raw_endpoints short_contig internal_window inclusive_floor original_floors output_label invalid_policy");
	}
	private static void coverageTiers(NcrnaFamily original,Path dir,Random random){
		final String[] names=Arrays.copyOf(original.modelNames,7);
		names[5]="mitochondria_coverage_fixture_first";names[6]="mitochondria_coverage_fixture_second";
		final byte[][] library=Arrays.copyOf(original.library,7);
		library[5]=randomBases(random,800);library[6]=randomBases(random,900);
		final String path=dir.resolve("coverage.fa").toString();final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{for(int i=0;i<library.length;i++){out.print(b.clear().append('>').append(names[i]).nl().append(library[i]).nl());}}
		finally{RrnaResourceIO.finish(out);}
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();options.resources.consensus=path;
		options.resources.lengthBounds=writeBounds(dir,"coverage.bounds.tsv",names,library,640,960,720,1080);
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("16S",options,new FamilyResourceResolver());
		check(loaded.library.length==7 && MitoSsuModelPolicy.coverage(loaded.name,loaded.modelNames),
			"Coverage selection may add a measured model count instead of requiring the old five additions");
		check(loaded.modelThresholds.minimum(5)==640 && loaded.modelThresholds.maximum(6)==1080
			&& loaded.modelThresholds.minimum(1)==0 && loaded.modelThresholds.maximum(1)==Integer.MAX_VALUE,
			"The unified loader must bind the added models' limits without imposing any on the original roster");
		for(int i=0;i<7;i++){check(loaded.rrnaEndpointFeatures.raw[i]==(i==4 || i>=5),
			"Coverage models keep raw ends while the original endpoint resource indices remain intact");}
		final byte[] gene=original.library[1];
		library[5]=new byte[gene.length+200];System.arraycopy(gene,0,library[5],0,gene.length);System.arraycopy(gene,0,library[5],gene.length,200);
		final LongHashSet seeds=seeds(gene);
		final TrnaKmerIndex index=new TrnaKmerIndex(library,9,false,0,0,0,1);
		check(index.shortlist(gene,1)[0]==5,"Fixture must put a passing added model ahead of the true original by K9 count");
		final NcrnaScavenger baseline=caller(original.library,original.modelNames,seeds),tiered=caller(library,names,seeds);
		final ArrayList<Orf> before=baseline.scavenge("original_first",gene,0,new ArrayList<int[]>());
		final ArrayList<Orf> after=tiered.scavenge("original_first",gene,0,new ArrayList<int[]>());
		check(before.size()==1 && after.size()==1 && before.get(0).toGff().equals(after.get(0).toGff())
			&& baseline.alignmentCount()==tiered.alignmentCount(),
			"A higher-ranked passing mito model must not change the original full call or its alignment count");
		final byte[] partial=Arrays.copyOfRange(gene,gene.length-1200,gene.length);
		library[5]=partial.clone();final LongHashSet partialSeeds=seeds(partial);
		final NcrnaScavenger partialBaseline=caller(original.library,original.modelNames,partialSeeds),partialTier=caller(library,names,partialSeeds);
		final ArrayList<Orf> oldPart=partialBaseline.scavenge("original_partial",partial,0,new ArrayList<int[]>());
		final ArrayList<Orf> newPart=partialTier.scavenge("original_partial",partial,0,new ArrayList<int[]>());
		check(oldPart.size()==1 && newPart.size()==1 && oldPart.get(0).ncrnaPartial!=0
			&& oldPart.get(0).toGff().equals(newPart.get(0).toGff()) && partialBaseline.alignmentCount()==partialTier.alignmentCount(),
			"Original physical-edge rescue must precede an exact full-length mitochondrial competitor");
		final byte[] mito=randomBases(random,800);library[5]=new byte[1600];
		System.arraycopy(mito,0,library[5],0,800);System.arraycopy(mito,0,library[5],800,800);library[6]=mito;
		final TrnaKmerIndex fallbackIndex=new TrnaKmerIndex(library,9,false,0,0,0,1);
		check(fallbackIndex.shortlist(mito,1)[0]==5,"Duplicated wrong-length model must lead the mito tier by shared K9 multiplicity");
		final NcrnaScavenger fallback=caller(library,names,seeds(mito));final ArrayList<Integer> attempted=new ArrayList<Integer>();
		final int[] finalRank={-1};
		fallback.setModelAttemptSink((n,s,p,a,z,rank,m,h,id,len,ok)->{attempted.add(m);if(m==6){finalRank[0]=rank;}});
		final ArrayList<Orf> rescued=fallback.scavenge("mito_fallback",mito,0,new ArrayList<int[]>());
		for(int m:attempted){check(m>=5,"An800-base edge window is below the original1200-base gate and must spend no original-model alignments");}
		check(rescued.size()==1 && rescued.get(0).trnaModel.equals(names[6]) && rescued.get(0).start==0 && rescued.get(0).stop==799
			&& attempted.indexOf(5)>=0 && attempted.indexOf(6)>attempted.indexOf(5) && finalRank[0]==2,
			"A failed leading mito model must fall through to the next ranked model without changing raw endpoints");
		lengthGates(original,dir,names,library,mito);
		System.out.println("MITO_ORIGINAL_FIRST_TIER_TEST_PASS dynamic_loader passing_competitor original_partial mito_fallback");
	}
	private static void lengthGates(NcrnaFamily original,Path dir,String[] names,byte[][] library,byte[] mito){
		final String fasta=dir.resolve("length_models.fa").toString();final ByteStreamWriter out=RrnaResourceIO.writer(fasta);final ByteBuilder b=new ByteBuilder();
		try{for(int i=0; i<library.length; i++){out.print(b.clear().append('>').append(names[i]).nl().append(library[i]).nl());}}
		finally{RrnaResourceIO.finish(out);}
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();options.resources.consensus=fasta;
		options.resources.lengthBounds=writeBounds(dir,"lengths.tsv",names,library,1280,1920,640,960);
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("16S",options,new FamilyResourceResolver());
		final NcrnaScavenger bounded=caller(library,names,seeds(mito));bounded.setModelThresholds(loaded.modelThresholds);
		final ArrayList<Integer> attempts=new ArrayList<Integer>();bounded.setModelAttemptSink((n,s,p,a,z,r,m,h,id,len,ok)->attempts.add(m));
		final ArrayList<Orf> calls=bounded.scavenge("bounded_fallback",mito,0,new ArrayList<int[]>());
		check(calls.size()==1 && calls.get(0).trnaModel.equals(names[6]) && !attempts.contains(5)
			&& bounded.mitoLengthSkips()>0 && bounded.mitoAlignmentCount()==1 && bounded.original16sAlignmentCount()==0,
			"An infeasible first model must spend zero alignments and permit the next model to recover the complete gene");
		final byte[] padded=new byte[3000];Arrays.fill(padded,(byte)'N');System.arraycopy(mito,0,padded,1000,800);
		// The duplicated fallback decoy can pass at.75 by extending800bp into N
		// padding (frozen replay9b2091dc). Isolate window geometry from that valid
		// first-accepted-model competition with an unrelated same-length decoy.
		final byte[][] paddedLibrary=library.clone();paddedLibrary[5]=randomBases(new Random(27221),1600);
		final NcrnaScavenger paddedCaller=caller(paddedLibrary,names,seeds(mito));paddedCaller.setModelThresholds(loaded.modelThresholds);
		final ArrayList<Orf> paddedCalls=paddedCaller.scavenge("padded_length",padded,0,new ArrayList<int[]>());
		check(paddedCalls.size()==1 && paddedCalls.get(0).start==1000 && paddedCalls.get(0).stop==1799 && paddedCalls.get(0).trnaModel.equals(names[6]),
			"A padded window exceeding max gene length must retain a valid internal gene span");
		final byte[] gene=original.library[1];final NcrnaScavenger baseline=caller(original.library,original.modelNames,seeds(gene));
		final NcrnaScavenger originalFirst=caller(library,names,seeds(gene));originalFirst.setModelThresholds(loaded.modelThresholds);
		final ArrayList<Orf> before=baseline.scavenge("bounded_original",gene,0,new ArrayList<int[]>()),after=originalFirst.scavenge("bounded_original",gene,0,new ArrayList<int[]>());
		check(before.size()==1 && after.size()==1 && before.get(0).toGff().equals(after.get(0).toGff())
			&& baseline.alignmentCount()==originalFirst.alignmentCount() && originalFirst.mitoAlignmentCount()==0,
			"Bound profiles must leave original-model calls and verification work unchanged");
		final String upper=writeBounds(dir,"upper.tsv",names,library,1280,1920,640,799);
		NcrnaModelThresholds.loadLengths(upper,loaded);
		final NcrnaScavenger overlong=caller(library,names,seeds(mito));overlong.setModelThresholds(loaded.modelThresholds);
		check(overlong.scavenge("overlong_gene",mito,0,new ArrayList<int[]>()).isEmpty() && overlong.mitoAlignmentCount()>0,
			"An exact800-base gene must fail a799-base upper limit after real verification");
		final String cutoff=dir.resolve("identity.tsv").toString();final ByteStreamWriter ids=RrnaResourceIO.writer(cutoff);
		try{ids.println("family\tmodel\tidpass\tidborderline");for(String name:names){ids.print(b.clear().append("16S\t").append(name).append("\t0.67\t0.66\n"));}}
		finally{RrnaResourceIO.finish(ids);}
		final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();families.add(loaded);NcrnaModelThresholds.load(cutoff,families);
		check(loaded.modelThresholds.pass(6)==.67f && loaded.modelThresholds.maximum(6)==799,
			"Later identity overrides must preserve already bound gene-length limits");
		final NcrnaModelThresholds prior=loaded.modelThresholds;
		final byte[][] wrong=library.clone();wrong[6]=Arrays.copyOf(mito,801);
		final String bad=writeBounds(dir,"wrong_length.tsv",names,wrong,1280,1920,640,960);
		boolean rejected=false;try{NcrnaModelThresholds.loadLengths(bad,loaded);}catch(IllegalArgumentException e){rejected=true;}
		check(rejected && loaded.modelThresholds==prior,"A changed consensus length must fail before changing live cutoffs");
		System.out.println("MITO_LENGTH_GATE_TEST_PASS loader prefilter fallback padded_window original_preserved upper_bound overrides atomic_failure");
	}
	private static String writeBounds(Path dir,String file,String[] names,byte[][] library,int lo5,int hi5,int lo6,int hi6){
		check(names.length==7 && library.length==7,"Length fixture binds exactly two added models");
		final String path=dir.resolve(file).toString();final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{out.println("model\tmodel_length\tmin_length\tmax_length");for(int i=5; i<7; i++){
			out.print(b.clear().append(names[i]).tab().append(library[i].length).tab().append(i==5 ? lo5 : lo6).tab().append(i==5 ? hi5 : hi6).nl());}}
		finally{RrnaResourceIO.finish(out);}return path;
	}
	private static void tierLimits(NcrnaFamily original,Path dir,Random random){
		final String[] properties={"bbtools.mito.maxModels","bbtools.mito.minHits","bbtools.mito.queryFraction","bbtools.mito.topFraction","bbtools.diagnostic.mitoTierOut"};
		final String[] saved=new String[properties.length];for(int i=0; i<saved.length; i++){saved[i]=System.getProperty(properties[i]);System.clearProperty(properties[i]);}
		try{
			final MitoSsuTierPolicy defaults=MitoSsuTierPolicy.configured();
			check(defaults.maxModels==3 && defaults.minHits==80 && defaults.queryFraction==0 && defaults.topFraction==0,
				"Cleared system properties must select the measured cap3/floor80/query0/top0 production defaults");
			final byte[] gene=randomBases(random,800),window=new byte[3000];Arrays.fill(window,(byte)'N');System.arraycopy(gene,0,window,1000,800);
			final byte[][] library=Arrays.copyOf(original.library,8);final String[] names=Arrays.copyOf(original.modelNames,8);
			for(int i=5; i<7; i++){library[i]=new byte[1600];System.arraycopy(gene,0,library[i],0,800);System.arraycopy(gene,0,library[i],800,800);names[i]="mitochondria_coverage_tier_"+i;}
			library[7]=gene;names[7]="mitochondria_coverage_tier_7";
			final NcrnaFamily family=new NcrnaFamily("16S",library,null,names,seeds(gene),17,1400,1800,9,1,false,0,0,0,1,0,1,.66f,.66f);family.reuseConsensusAlignment=true;
			final String bounds=dir.resolve("tier_bounds.tsv").toString();final ByteStreamWriter out=RrnaResourceIO.writer(bounds);final ByteBuilder b=new ByteBuilder();
			try{out.println("model\tmodel_length\tmin_length\tmax_length");for(int i=5; i<8; i++){
				out.print(b.clear().append(names[i]).tab().append(library[i].length).tab().append(i==7?640:1601).tab().append(i==7?960:2000).nl());}}
			finally{RrnaResourceIO.finish(out);}NcrnaModelThresholds.loadLengths(bounds,family);
			final NcrnaScavenger bare=caller(library,names,seeds(gene));bare.setModelThresholds(family.modelThresholds);
			final ArrayList<Orf> bareCalls=bare.scavenge("bare_default_third",window,0,new ArrayList<int[]>());
			check(bareCalls.size()==1 && bareCalls.get(0).trnaModel.equals(names[7]) && bareCalls.get(0).start==1000 && bareCalls.get(0).stop==1799,
				"The actual caller with no tier properties must reach the third model after two failed decoys and retain raw endpoints");
			System.out.println("MITO_TIER_DEFAULT_TEST_PASS cap3 floor80 query0 top0 actual_third_model_success");
			final String log=dir.resolve("tier_log.tsv").toString();System.setProperty(properties[4],log);
			System.setProperty(properties[0],"2147483647");
			final NcrnaScavenger replay=caller(library,names,seeds(gene));replay.setModelThresholds(family.modelThresholds);
			final ArrayList<Orf> accepted=replay.scavenge("uncapped",window,0,new ArrayList<int[]>());
			check(accepted.size()==1 && accepted.get(0).trnaModel.equals(names[7]),"Fixture must require the third model after two length-ineligible alignment results");
			System.setProperty(properties[0],"2");
			final NcrnaScavenger capped=caller(library,names,seeds(gene));capped.setModelThresholds(family.modelThresholds);
			check(capped.scavenge("capped",window,0,new ArrayList<int[]>()).isEmpty() && capped.mitoAlignmentCount()==4,
				"Two distinct models retain their full/edge verifier calls but cannot fall through to a third model");
			System.setProperty(properties[1],"65536");
			final NcrnaScavenger filtered=caller(library,names,seeds(gene));filtered.setModelThresholds(family.modelThresholds);
			check(filtered.scavenge("entry_rejected",window,0,new ArrayList<int[]>()).isEmpty() && filtered.mitoAlignmentCount()==0,
				"An entry-floor failure must spend zero mitochondrial alignments");
			System.setProperty(properties[1],"1000");
			final NcrnaScavenger successor=caller(library,names,seeds(gene));successor.setModelThresholds(family.modelThresholds);
			check(successor.scavenge("entry_successor",gene,0,new ArrayList<int[]>()).isEmpty() && successor.mitoAlignmentCount()==0 && successor.mitoLengthSkips()>0,
				"An infeasible high-count model must not admit a feasible successor below the entry floor");
			System.setProperty(properties[0],"3");System.setProperty(properties[1],"0");System.setProperty(properties[3],"1");
			final NcrnaScavenger relative=caller(library,names,seeds(gene));relative.setModelThresholds(family.modelThresholds);
			check(relative.scavenge("relative_rejected",window,0,new ArrayList<int[]>()).isEmpty() && relative.mitoAlignmentCount()==4,
				"Relative-best filtering must reject the lower-count true model after two actual decoy verifications");
			final TrnaKmerIndex boundaryIndex=new TrnaKmerIndex(library,9,false,0,0,0,1);boundaryIndex.shortlist(window,8);
			final int topCount=boundaryIndex.lastSharedCount(5),trueCount=boundaryIndex.lastSharedCount(7);
			float inclusiveFraction=(float)trueCount/topCount;
			while(Math.ceil((double)inclusiveFraction*topCount)>trueCount){inclusiveFraction=Math.nextDown(inclusiveFraction);}
			check(new MitoSsuTierPolicy(3,0,0,inclusiveFraction).subsequentFloor(topCount,boundaryIndex.lastUniqueQueryKmers())==trueCount,
				"The relative fixture must put the true model exactly on the integral shared-count cutoff despite float rounding");
			System.setProperty(properties[3],Float.toString(inclusiveFraction));
			final NcrnaScavenger boundary=caller(library,names,seeds(gene));boundary.setModelThresholds(family.modelThresholds);
			check(boundary.scavenge("relative_boundary",window,0,new ArrayList<int[]>()).size()==1,
				"A model exactly at the inclusive relative-score boundary must remain eligible in the actual caller");
			final byte[] background=randomBases(new Random(31255),window.length);System.arraycopy(gene,0,background,1000,gene.length);
			System.setProperty(properties[3],"0");System.setProperty(properties[2],"1");
			final NcrnaScavenger queryRejected=caller(library,names,seeds(gene));queryRejected.setModelThresholds(family.modelThresholds);
			check(queryRejected.scavenge("query_rejected",background,0,new ArrayList<int[]>()).isEmpty() && queryRejected.mitoAlignmentCount()==0,
				"A query-unique fraction can reject entry independently of the absolute floor");
			System.setProperty(properties[2],"0.25");
			final NcrnaScavenger queryAdmitted=caller(library,names,seeds(gene));queryAdmitted.setModelThresholds(family.modelThresholds);
			queryAdmitted.scavenge("query_admitted",background,0,new ArrayList<int[]>());
			check(queryAdmitted.mitoAlignmentCount()>0,"Lowering only the query fraction must admit actual verification work on the same window");
			MitoSsuTierLog.closeConfigured();
			final fileIO.ByteFile in=RrnaResourceIO.open(log);final parse.LineParser1 p=new parse.LineParser1('\t');int rows=0;
			try{RrnaResourceIO.header(in,MitoSsuTierLog.HEADER,log);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);rows++;
				check(p.terms()==18 && p.parseInt(5)>0,"Observer must retain unique query9-mer counts and its complete schema");
				if(p.termEquals("uncapped",0)){check(p.parseInt(9)==3 && p.parseInt(14)==3 && p.termEquals("ACCEPT",16),"Uncapped trace must expose third-attempt acceptance for calibration");}
				else if(p.termEquals("capped",0)){check(p.parseInt(9)==2 && p.parseInt(10)==4 && p.parseInt(17)==4 && p.termEquals("REJECT_CAP",16),"Cap trace must separate distinct models from aligner calls");}
				else if(p.termEquals("relative_rejected",0)){check(p.parseInt(9)==2 && p.termEquals("REJECT_RELATIVE",16),"Relative rejection must be distinct from the model cap");}
				else if(p.termEquals("relative_boundary",0)){check(p.termEquals("ACCEPT",16),"Inclusive relative eligibility must retain its observed acceptance");}
				else if(p.termEquals("query_admitted",0)){check(p.parseInt(9)>0,"Query threshold admission must appear in the observation");}
				else{check((p.termEquals("entry_rejected",0) || p.termEquals("entry_successor",0) || p.termEquals("query_rejected",0))
					&& p.parseInt(9)==0 && p.termEquals("REJECT_ENTRY",16),"Entry-rejected windows remain observable even with zero attempts");}
			}}finally{check(!in.close(),"Tier observer fixture read must finish");}
			check(rows==8,"Exactly one observer row is required per window reaching the mitochondrial tier");
			final MitoSsuTierPolicy policy=new MitoSsuTierPolicy(2,10,.25f,.75f);
			check(policy.entryFloor(100)==25 && policy.subsequentFloor(100,100)==75,"Query-based entry and relative-best subsequent thresholds have distinct meanings");
		}finally{MitoSsuTierLog.closeConfigured();for(int i=0; i<saved.length; i++){if(saved[i]==null){System.clearProperty(properties[i]);}else{System.setProperty(properties[i],saved[i]);}}}
		System.out.println("MITO_TIER_LIMIT_TEST_PASS uncapped_observer cap_distinct_models preserved_edge_attempts zero_alignment_entry_gate");
	}
	private static LongHashSet seeds(byte[] gene){
		check(gene.length>400,"Fixture needs two separated complete forward seed keys");
		final LongHashSet result=new LongHashSet(8);result.add(key(gene,80));result.add(key(gene,350));return result;
	}
	private static NcrnaScavenger caller(final byte[][] library, final String[] names, final LongHashSet seeds){
		final NcrnaScavenger result=new NcrnaScavenger(library,null,names,seeds,17,1400,1800,9,1,false,0,0,0,1,0,1,.66f,.66f);
		result.family="16S";result.outputType=ProkObject.r16S;result.reuseConsensusAlignment=true;result.trimAlignmentExtent=false;
		result.prok16sPartialMinLen=1200;result.physicalEdgeMinLen=1200;result.minKmerHits=2;result.scavengePass2=false;
		return result;
	}
	private static byte[] randomBases(final Random random, final int length){
		assert(length>0):"Fixture consensus must have a positive complete length";
		final byte[] bases=new byte[length],alphabet={'A','C','G','T'};for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}return bases;
	}
	private static long key(final byte[] bases, final int start){
		assert(start>=0 && start+17<=bases.length):"Fixture seeds must be genuine full forward17mers";
		long word=0;for(int i=start; i<start+17; i++){word=(word<<2)|RrnaResourceIO.base(bases[i]);}return word;
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
