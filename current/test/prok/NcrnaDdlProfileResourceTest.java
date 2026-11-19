package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Locale;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import dna.AminoAcid;
import kmer.AbstractKmerTable;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlInterleave;
import prok.NcrnaDdlDualKeys;
import prok.NcrnaDdlPartitionKeys;
import prok.NcrnaDdlProfile;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaDdlQueryKeys;
import prok.NcrnaDdlOrf;
import prok.NcrnaDdlShortlist;
import prok.NcrnaDdlRanks;
import prok.NcrnaAlignmentWork;
import prok.NcrnaCpuProfile;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.TrnaKmerIndex;
import prok.GeneCaller;
import prok.GeneModel;
import prok.NcrnaFamily;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Small resource round-trip, score/order parity and explicit binding faults.
 * No caller, genome or biological accuracy measurement.
 * @author Brian, Raiden
 */
public final class NcrnaDdlProfileResourceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlProfileResourceTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String verify=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("verify")){
					check(b!=null && !b.isEmpty(), "verify requires the actual tiny builder output path");
					verify=b;
				}
				else{check(parser.parse(arg, a, b), "Unknown resource fixture option: "+arg);}
			}
			check((verify==null ? parser.out1!=null : parser.out1==null) && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use either out=fresh_fixture_directory or verify=tiny_builder_output.tsv without input, overwrite or sampling");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Resource and rank fixtures require assertions");
			DynamicDemiLog.setExponent(4);
			if(verify!=null){
				verifyBuilder(verify);
				pp.outstream.println("NCRNA_DDL_PROFILE_CLI_ROUNDTRIP_PASS");
			}else{
				final File root=new File(parser.out1).getAbsoluteFile();
				check(root.mkdir(), "Fixture directory must be fresh");
				run(root.getPath());
				pp.outstream.println("NCRNA_DDL_PROFILE_RESOURCE_TEST_PASS");
			}
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(final String root){
		inverted();
		final byte[][] sequences={bytes("aaaaaaaaaaaaaaaaaaaa"), bytes("CCCCCCCCCCCCCCCCCCCC")};
		final String[] names={"cobalamin_RF00174_α n=2", "cobalamin_RF00174_beta n=1"};
		final NcrnaDdlProfile[] profiles=profiles(sequences, 2);
		final String path=root+"/good.tsv";
		NcrnaDdlProfileResource.write(path, profiles, names, sequences);
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, sequences);
		check(resource.models()==2 && resource.members(0)==2 && resource.members(1)==1,
			"All model rows retain their independent voter counts");
		for(int model=0; model<profiles.length; model++){
			for(int bucket=0; bucket<NcrnaDdlProfileResource.BUCKETS; bucket++){
				if(resource.key(model, bucket)!=profiles[model].keys()[bucket]){
					throw new AssertionError("Resource changed an absolute key: model="+model+", bucket="+bucket);
				}
			}
		}
		parity(resource, profiles);
		reject(()->NcrnaDdlProfileResource.write(path, profiles, names, sequences), "Output exists");
		final byte[][] changed={sequences[0].clone(), sequences[1].clone()};
		changed[0][0]='G';
		reject(()->NcrnaDdlProfileResource.load(path, names, changed), "consensus content");
		reject(()->NcrnaDdlProfileResource.load(path, new String[]{names[1], names[0]}, sequences), "model ID");
		reject(()->NcrnaDdlProfileResource.load(path, new String[]{names[0], names[0]}, sequences), "unique");
		final String[] faults={"header", "missing", "extra", "order", "name", "sequence", "members", "negative", "overflow", "short", "duplicate_id"};
		final String[] reasons={"header", "Missing profile", "beyond", "original model order", "model ID", "consensus content",
			"member votes", "unsigned16bits", "unsigned16bits", "complete bucket vector", "model ID"};
		check(faults.length==reasons.length, "Each malformed resource needs its intended diagnostic");
		for(int i=0; i<faults.length; i++){
			final String bad=root+'/'+faults[i]+".tsv";
			fault(bad, faults[i], profiles, names, sequences);
			final String reason=reasons[i];
			reject(()->NcrnaDdlProfileResource.load(bad, names, sequences), reason);
		}
		final char before=resource.key(0, 0);
		profiles[0].keys()[0]++;
		check(resource.key(0, 0)==before, "Loaded keys must not alias the producer's mutable profile arrays");
		literalKeys(root, names, sequences);
		callerIntegration(resource, names);
		physicalEdgeIntegration(root);
		dual(root, names, sequences);
		compact(root, names, sequences);
		if(NcrnaCpuProfile.enabled()){profiledWindows(root);}
		globalIndependence(root, resource);
	}

	/** The finite recipe enables the real profiler before class initialization.
	 * A seeded early reject is one zero-work window; a physical partial includes
	 * both failed primary work and its retry inside one scheduled-window unit. */
	private static void profiledWindows(final String root){
		final byte[][] small={bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("CCCCCCCCCCCCCCCCCCCC")};
		final String[] smallNames={"cobalamin_RF00174_zero0", "cobalamin_RF00174_zero1"};
		final NcrnaScavenger zero=scavenger(small, smallNames);
		zero.minKmerHits=Integer.MAX_VALUE;
		final NcrnaCpuProfile zeroProfile=NcrnaCpuProfile.create(1);
		final ArrayList<Orf> empty=observed(zeroProfile, zero, small[0]);
		final NcrnaAlignmentWork zeroWork=zeroProfile.alignmentWork(0);
		check(empty.isEmpty() && zero.windowCount()==1 && zeroWork.windows==1 && zeroWork.logical==0
			&& zeroWork.physical==0 && zeroWork.windowCount(true, 0)==1 && zeroWork.windowCount(false, 0)==1,
			"Every scheduled seeded rejection must contribute one zero bin, not disappear from the workload denominator");
		writeProfile(zeroProfile, zero, root+"/zero_cpu.tsv");
		final byte[][] models={Prok16sPartialRescueTest.random(1000), Prok16sPartialRescueTest.random(1000, 731)};
		final String[] names={"cobalamin_RF00174_work0", "cobalamin_RF00174_work1"};
		final byte[] fragment=Arrays.copyOfRange(models[0], 200, 850);
		final String legacyPath=root+"/work_legacy.tsv", dualPath=root+"/work_dual.tsv";
		NcrnaDdlProfileResource.write(legacyPath, profiles(models, 1), names, models);
		NcrnaDdlProfileResource.write(dualPath, dualProfiles(models), names, models, true);
		final String compact512=root+"/work_compact512.tsv", compact1024=root+"/work_compact1024.tsv";
		NcrnaDdlProfileResource.writeCompact(compact512, compactProfiles(models, 512), names, models, 512);
		NcrnaDdlProfileResource.writeCompact(compact1024, compactProfiles(models, 1024), names, models, 1024);
		final NcrnaDdlProfileResource[] resources={null,
			NcrnaDdlProfileResource.load(legacyPath, names, models), NcrnaDdlProfileResource.load(dualPath, names, models),
			NcrnaDdlProfileResource.load(compact512, names, models), NcrnaDdlProfileResource.load(compact1024, names, models)};
		for(int mode=0; mode<resources.length; mode++){
			final NcrnaDdlProfileResource resource=resources[mode];
			final NcrnaScavenger plain=edgeScavenger(models, names, resource==null || !resource.dual());
			final NcrnaScavenger worker=edgeScavenger(models, names, resource==null || !resource.dual());
			if(resource!=null){
				plain.useDdlProfiles(resource, 50);
				worker.useDdlProfiles(resource, 50);
			}
			final ArrayList<Orf> expected=plain.scavenge("observed", fragment, 0, new ArrayList<int[]>());
			final NcrnaCpuProfile profile=NcrnaCpuProfile.create(1);
			final ArrayList<Orf> calls=observed(profile, worker, fragment);
			final NcrnaAlignmentWork work=profile.alignmentWork(0);
			check(calls.size()==1 && expected.size()==1 && calls.get(0).ncrnaPartial==3
				&& calls.get(0).toGff().equals(expected.get(0).toGff()),
				"Observing the complete ordinary/hybrid/dual edge window must preserve the exact synthetic call and metadata");
			check(worker.windowCount()==1 && work.windows==1 && work.logical>=2 && work.logical==work.physical
				&& work.windowLogical==work.logical && work.windowPhysical==work.physical,
				"Failed primary and physical-edge retry kernels must remain in one window with no unassigned work");
			writeProfile(profile, worker, root+"/edge_cpu_"+mode+".tsv");
		}
		System.err.println("NCRNA_DDL_WINDOW_WORK_TEST_PASS");
	}

	private static ArrayList<Orf> observed(final NcrnaCpuProfile profile, final NcrnaScavenger worker, final byte[] sequence){
		check(profile!=null, "Real window accounting fixture requires the profiling JVM property before class initialization");
		profile.begin(0, NcrnaCpuProfile.SCAVENGE);
		try{return worker.scavenge("observed", sequence, 0, new ArrayList<int[]>());}
		finally{profile.end();}
	}

	private static void writeProfile(final NcrnaCpuProfile profile, final NcrnaScavenger worker, final String path){
		check(profile!=null && worker!=null, "A profile report needs the actual observed worker and its completed counters");
		profile.activity(0, worker.alignmentCount(), worker.alignedBases(), worker.hbmScoreCalls(), worker.hbmBasesScored(), worker.kmerHitCount(), worker.windowCount());
		profile.verification(0, worker.alignmentCpuNanos(), worker.shortlistCpuNanos(), worker.hbmCpuNanos(), worker.modelLengthChecks(), worker.modelLengthSkips(), worker.modelSpanRejects());
		profile.write(path, new String[]{"cobalamin_rf00174"});
	}

	/** Independent single-length queries check the serialized dual half scores.
	 * Complete zero-score order, native-index absence and ordinary-worker reuse
	 * are distinct from the hybrid's positive-only proposal union. */
	private static void dual(final String root, final String[] names, final byte[][] sequences){
		final NcrnaDdlProfile[] profiles=dualProfiles(sequences);
		final String path=root+"/dual.tsv";
		NcrnaDdlProfileResource.write(path, profiles, names, sequences, true);
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, sequences);
		check(resource.dual() && resource.buckets()==4096 && resource.models()==2,
			"The explicit V2 header must retain the selected partitioned geometry, not inherit V1 width");
		final char[][] shorts=new char[2][], longs=new char[2][];
		for(int model=0; model<2; model++){
			shorts[model]=Arrays.copyOfRange(profiles[model].keys(), 0, 2048);
			longs[model]=Arrays.copyOfRange(profiles[model].keys(), 2048, 4096);
			for(int bucket=0; bucket<4096; bucket++){
				check(resource.key(model, bucket)==profiles[model].keys()[bucket], "V2 round-trip must retain every half/bucket key");
			}
		}
		final NcrnaDdlQueryKeys shortQuery=new NcrnaDdlQueryKeys(2048, 9), longQuery=new NcrnaDdlQueryKeys(2048, 10);
		final NcrnaDdlProfileResource.Scorer scorer=resource.newScorer();
		final NcrnaDdlShortlist order=new NcrnaDdlShortlist(resource, 384);
		final float[] scores=new float[2];
		int[] buffer=null;
		for(final String text:new String[]{"AAAAAAAAAAAAAAAAAAAA", "cccccccccccccccccccc", "", "NNNNNNNNNNNN", "AAAAACAAAA"}){
			final byte[] bases=bytes(text);
			shortQuery.build(bases, null, false, false);
			longQuery.build(bases, null, false, true);
			scorer.score(bases, scores);
			for(int model=0; model<2; model++){
				final float expected=.5f*shortQuery.equalBuckets(shorts[model])+longQuery.equalBuckets(longs[model]);
				check(scores[model]==expected, "One dual posting traversal must equal independent short/long scalar bucket queries");
			}
			final int[] observed=order.build(bases, null);
			check(buffer==null || buffer==observed, "Sketch-only queries must reuse their borrowed model buffer");
			buffer=observed;
			final int first=scores[1]>scores[0] ? 1 : 0;
			check(order.size()==2 && observed[0]==first && observed[1]==1-first
				&& order.nativeCount(first)==-1 && order.sketchFirst(first),
				"Sketch-only order includes zero-score models, resolves ties by ordinal and leaves native support unmeasured");
		}
		final NcrnaDdlInterleave many=new NcrnaDdlInterleave(508);
		final float[] values=new float[508];
		values[7]=1.5f;
		values[40]=1.5f;
		for(final int cap:new int[]{50,384}){
			check(many.buildSketchCaller(values, cap)==cap && many.modelAt(0)==7 && many.modelAt(1)==40
				&& many.modelAt(2)==0 && many.modelAt(cap-1)==cap-1,
				"Both real caps must apply to the complete ordinal-tied sketch order, including its zero-score tail");
		}
		Arrays.fill(values, 0);
		check(many.buildSketchCaller(values, 50)==50 && many.modelAt(0)==0 && many.modelAt(49)==49,
			"A zero-score replacement keeps the full saved order and clears previous model ranks");
		final byte[][] callerSequences={bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("CCCCCCCCCCCCCCCCCCCC")};
		final String baseline=scavenger(callerSequences, names).scavenge("fixture", bytes("AAAAAAAAAAAAAAAAAAAA"), 0, new ArrayList<int[]>()).get(0).toGff();
		runOwnership(resource, callerSequences, names, bytes("AAAAAAAAAAAAAAAAAAAA"), baseline);
		dualEdge(root);
	}

	private static NcrnaDdlProfile[] dualProfiles(final byte[][] sequences){
		check(sequences.length==2, "The dual integration fixture binds two original models");
		final NcrnaDdlDualKeys keys=new NcrnaDdlDualKeys(9, 10, 4096, true);
		final NcrnaDdlProfile[] profiles=new NcrnaDdlProfile[2];
		for(int model=0; model<2; model++){
			keys.build(sequences[model]);
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(4096);
			votes.add(keys.memberKeys());
			profiles[model]=votes.finish();
		}
		return profiles;
	}

	/** Compare the serialized caller scorer to direct raw bucket/key membership,
	 * not to another posting traversal. Check the complete saved-rank tie rule. */
	private static void compactParity(final NcrnaDdlProfileResource resource, final NcrnaDdlProfile[] profiles){
		check(resource.compact() && resource.dual() && profiles.length==2,
			"Both compact widths must select sketch-only ordering for the bound two-model fixture");
		final NcrnaDdlPartitionKeys query=new NcrnaDdlPartitionKeys(7, resource.buckets());
		final NcrnaDdlProfileResource.Scorer first=resource.newScorer(), second=resource.newScorer();
		final NcrnaDdlShortlist order=new NcrnaDdlShortlist(resource, 384);
		final float[] scores=new float[2], other=new float[2], expected=new float[2];
		final int[] rank=new int[4];
		final boolean[] acceptable=new boolean[2];
		final java.util.Random random=new java.util.Random(20261114);
		int[] buffer=null;
		for(int sample=0; sample<64; sample++){
			final byte[] bases=sample==0 ? bytes("AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCC")
				: sample==1 ? bytes("aaaaaaaNNNNNNNNNNRcccccccccc") : sample==2 ? bytes("") : new byte[sample*7];
			if(sample>2){
				for(int i=0; i<bases.length; i++){bases[i]=(byte)"ACGTN".charAt(random.nextInt(5));}
			}
			query.build(bases);
			Arrays.fill(expected, 0);
			for(int i=0; i<query.keyCount(); i++){
				final int packed=query.packedKey(i), bucket=packed>>>16;
				for(int model=0; model<2; model++){
					if(profiles[model].keys()[bucket]==(char)packed){expected[model]++;}
				}
			}
			second.score(bytes("CCCCCCCCCCCCCCCCCCCC"), other);
			first.score(bases, scores);
			check(Arrays.equals(scores, expected),
				"Compact posting scores must equal direct unique raw-key membership, including empty/undefined/repeated words");
			final int[] observed=order.build(bases, null);
			check(order.size()==2 && (buffer==null || observed==buffer),
				"Compact ordering retains zero-score models and reuses the same borrowed output buffer");
			buffer=observed;
			for(int model=0; model<2; model++){
				acceptable[model]=true;
				NcrnaDdlRanks.best(expected, acceptable, rank);
				check(observed[rank[0]-1]==model && order.nativeCount(model)==-1,
					"Caller order must equal the saved descending-score/original-ordinal rank, with native support unmeasured");
				acceptable[model]=false;
			}
		}
	}

	private static NcrnaDdlProfile[] compactProfiles(final byte[][] sequences, final int buckets){
		assert(sequences.length==2) : "The resource fixture has two independent literal member/consensus pairs";
		final NcrnaDdlPartitionKeys codec=new NcrnaDdlPartitionKeys(7, buckets);
		final NcrnaDdlProfile[] result=new NcrnaDdlProfile[2];
		for(int model=0; model<2; model++){
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(buckets);
			votes.add(codec.memberKeys(sequences[model]));
			result[model]=votes.finishPartition(codec.consensusKeys(sequences[model]));
		}
		return result;
	}

	private static void compact(final String root, final String[] names, final byte[][] sequences){
		final byte[][] models={Prok16sPartialRescueTest.random(1000), Prok16sPartialRescueTest.random(1000, 731)};
		final byte[] fragment=Arrays.copyOfRange(models[0], 200, 850);
		// Resource binding accepts lowercase consensus bytes; the ordinary caller
		// fixture uses normalized models, as the production library loader does.
		final byte[][] callerSequences={bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("CCCCCCCCCCCCCCCCCCCC")};
		final ArrayList<Orf> baselineCalls=scavenger(callerSequences, names).scavenge("fixture", callerSequences[0], 0, new ArrayList<int[]>());
		check(baselineCalls.size()==1, "The normalized A20/C20 ordinary caller fixture must retain its one exact A20 baseline");
		final String baseline=baselineCalls.get(0).toGff();
		for(final int buckets:new int[]{512,1024}){
			final NcrnaDdlProfile[] profiles=compactProfiles(sequences, buckets);
			final String path=root+"/compact_"+buckets+".tsv";
			NcrnaDdlProfileResource.writeCompact(path, profiles, names, sequences, buckets);
			final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, sequences);
			check(resource.buckets()==buckets && resource.members(0)==1 && resource.members(1)==1,
				"The compact header must preserve its total width and original voter counts");
			for(int model=0; model<2; model++){
				for(int bucket=0; bucket<buckets; bucket++){
					check(resource.key(model, bucket)==profiles[model].keys()[bucket],
						"Compact serialization must preserve every member and shared-consensus bucket");
				}
			}
			compactParity(resource, profiles);
			runOwnership(resource, callerSequences, names, callerSequences[0], baseline);
			final String edgePath=root+"/compact_edge_"+buckets+".tsv";
			NcrnaDdlProfileResource.writeCompact(edgePath, compactProfiles(models, buckets), names, models, buckets);
			final NcrnaScavenger caller=edgeScavenger(models, names, false);
			caller.useDdlProfiles(NcrnaDdlProfileResource.load(edgePath, names, models), 384);
			final ArrayList<Orf> calls=caller.scavenge("physical", fragment, 0, new ArrayList<int[]>());
			final ArrayList<Orf> ordinary=edgeScavenger(models, names).scavenge("physical", fragment, 0, new ArrayList<int[]>());
			check(!caller.hasNativeRankingIndex() && caller.modelPatience()==8 && calls.size()==1 && ordinary.size()==1
				&& calls.get(0).ncrnaPartial==3
				&& calls.get(0).toGff().equals(ordinary.get(0).toGff()+";ddl_stream=sketch_only;ddl_native_shared=NA"),
				"Compact scoring must preserve the ordinary physical-edge verifier/patience and allocate no native index");
		}
		reject(()->NcrnaDdlProfileResource.writeCompact(root+"/wrong_width.tsv", compactProfiles(sequences, 512), names, sequences, 1024), "declared width");
		reject(()->NcrnaDdlProfileResource.writeCompact(root+"/unsupported.tsv", compactProfiles(sequences, 512), names, sequences, 2048), "512or1024");
		System.err.println("NCRNA_DDL_COMPACT_RESOURCE_TEST_PASS");
	}

	private static void dualEdge(final String root){
		final byte[][] models={Prok16sPartialRescueTest.random(1000), Prok16sPartialRescueTest.random(1000, 731)};
		final String[] names={"cobalamin_RF00174_edge0", "cobalamin_RF00174_edge1"};
		final String path=root+"/dual_physical.tsv";
		NcrnaDdlProfileResource.write(path, dualProfiles(models), names, models, true);
		final NcrnaDdlProfileResource profiles=NcrnaDdlProfileResource.load(path, names, models);
		final byte[] fragment=Arrays.copyOfRange(models[0], 200, 850);
		final String baseline=edgeScavenger(models, names).scavenge("physical", fragment, 0, new ArrayList<int[]>()).get(0).toGff();
		for(final int cap:new int[]{50,384}){
			final NcrnaScavenger caller=edgeScavenger(models, names, false);
			caller.useDdlProfiles(profiles, cap);
			check(!caller.hasNativeRankingIndex(), "Sketch-only edge retries must not allocate a native index");
			final ArrayList<Orf> calls=caller.scavenge("physical", fragment, 0, new ArrayList<int[]>());
			check(calls.size()==1 && calls.get(0).ncrnaPartial==3
				&& calls.get(0).toGff().equals(baseline+";ddl_stream=sketch_only;ddl_native_shared=NA"),
				"Sketch-only caps preserve physical-edge verification and report unmeasured native support honestly");
		}
	}

	private static void callerIntegration(final NcrnaDdlProfileResource resource, final String[] names){
		final byte[][] sequences={bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("CCCCCCCCCCCCCCCCCCCC")};
		final TrnaKmerIndex index=new TrnaKmerIndex(sequences, 9, true, 2f, 0f, 0f, 0);
		final NcrnaDdlShortlist order=new NcrnaDdlShortlist(resource, 50);
		final byte[] query=bytes("AAAAAAAAAAAAAAAAAAAA");
		final int[] buffer=order.build(query, index);
		check(order.size()==1 && buffer[0]==0 && order.nativeCount(0)==12 && !order.sketchFirst(0),
			"The literal A20 query has twelve native9mer occurrences for model0 and no eligible model1");
		check(order.build(bytes(""), index)==buffer && order.size()==0,
			"An empty replacement must expose zero live models while retaining the same borrowed array");
		check(order.build(bytes("CCCCCCCCCCCCCCCCCCCC"), index)==buffer && order.size()==1 && buffer[0]==1,
			"The next query replaces the live prefix without allocating another model array");
		check(order.build(bytes("AAAAACAAAA"), index)==buffer && order.size()==1 && buffer[0]==0
			&& order.nativeCount(0)==0 && order.sketchFirst(0),
			"A central mismatch destroys both native9mers but the forward masked10mer must still propose model0");
		final ArrayList<Orf> ordinary=scavenger(sequences, names).scavenge("fixture", query, 0, new ArrayList<int[]>());
		check(ordinary.size()==1 && !(ordinary.get(0) instanceof NcrnaDdlOrf),
			"Disabled hybrid caller must retain its ordinary candidate class");
		final String baseline=ordinary.get(0).toGff();
		for(final int cap:new int[]{50,384}){
			final NcrnaScavenger caller=scavenger(sequences, names);
			caller.useDdlProfiles(resource, cap);
			final ArrayList<Orf> calls=caller.scavenge("fixture", query, 0, new ArrayList<int[]>());
			check(calls.size()==1 && calls.get(0) instanceof NcrnaDdlOrf,
				"Both experimental caps must use the same ordinary verifier and retain the one literal call");
			final ArrayList<Orf> selected=NcrnaScavenger.selectBestNonOverlapping(calls);
			final String observed=selected.get(0).toGff();
			check(observed.equals(baseline+";ddl_stream=native;ddl_native_shared=12"),
				"Proposal metadata must survive selection while preserving every ordinary GFF field");
			final long before=caller.alignmentCount();
			check(caller.scavenge("empty_union", bytes("GGGGGGGGGGGGGGGGGGGG"), 0, new ArrayList<int[]>()).isEmpty()
				&& caller.alignmentCount()==before,
				"A seeded window with neither proposal stream must align no stale entries from the previous live prefix");
		}
		check(scavenger(sequences, names).scavenge("fixture", query, 0, new ArrayList<int[]>()).get(0).toGff().equals(baseline),
			"An ordinary caller after experimental workers must retain exact disabled-path output");
		runOwnership(resource, sequences, names, query, baseline);
	}

	/** Reuse one cached family across ordinary, experimental and ordinary workers.
	 * Installing on one run must never mutate the shared descriptor for the next. */
	private static void runOwnership(final NcrnaDdlProfileResource resource, final byte[][] sequences,
			final String[] names, final byte[] query, final String baseline){
		final ArrayList<NcrnaFamily> previous=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			final LongHashSet seeds=new LongHashSet(16);
			seeds.add(0);
			final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", sequences, null, names, seeds,
				17, 18, 20, 9, 384, true, 2f, 0f, 0f, 0, 0f, 20f, 1f, 1f);
			family.modelPatience=8;
			family.physicalEdgeMinLen=18;
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(family);
			for(int run=0; run<3; run++){
				final GeneCaller caller=new GeneCaller(60, 10, 10, -1f, -1f, -1f, -1f, -1f, new GeneModel(false));
				if(run==1){caller.setCobalaminDdl(resource, 50);}
				final NcrnaScavenger worker=caller.ncrnaScavenger(0);
				check(worker.hasNativeRankingIndex()==(run!=1 || !resource.dual()) && worker.modelPatience()==8,
					"Only sketch-only RF00174 skips native-index construction; every worker keeps shipped patience");
				final ArrayList<Orf> calls=worker.scavenge("fixture", query, 0, new ArrayList<int[]>());
				check(calls.size()==1 && (calls.get(0) instanceof NcrnaDdlOrf)==(run==1),
					"A later caller sharing the same cached family must not inherit another run's experiment");
				final String metadata=resource.dual() ? ";ddl_stream=sketch_only;ddl_native_shared=NA" : ";ddl_stream=native;ddl_native_shared=12";
				check(calls.get(0).toGff().equals(run==1 ? baseline+metadata : baseline),
					"Repeated caller instances must preserve ordinary output around the explicit experimental run");
			}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.addAll(previous);
		}
	}

	private static NcrnaScavenger scavenger(final byte[][] sequences, final String[] names){
		check(sequences.length==2 && names.length==2, "The integration fixture needs the bound A20/C20 model roster");
		final LongHashSet seeds=new LongHashSet(16);
		seeds.add(0);
		seeds.add(0x155555555L);
		seeds.add(0x2aaaaaaaaL);//G17 under the shared seed index's A/C/G/T encoding.
		final NcrnaScavenger caller=new NcrnaScavenger(sequences, null, names, seeds,
			17, 18, 20, 9, 384, true, 2f, 0f, 0f, 0, 0f, 20f, 1f, 1f);
		caller.family="cobalamin_rf00174";
		caller.setShortlistPolicy(384, true, 2f, 0f, 0f, 0, 8);
		caller.physicalEdgeMinLen=18;
		return caller;
	}

	/** Production edge policy keeps both original model ends in the clipped route. */
	private static void physicalEdgeIntegration(final String root){
		final byte[][] models={Prok16sPartialRescueTest.random(1000), Prok16sPartialRescueTest.random(1000, 731)};
		final String[] names={"cobalamin_RF00174_edge0", "cobalamin_RF00174_edge1"};
		final String path=root+"/physical_profiles.tsv";
		NcrnaDdlProfileResource.write(path, profiles(models, 1), names, models);
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, models);
		final byte[] fragment=Arrays.copyOfRange(models[0], 200, 850);
		final NcrnaScavenger primary=edgeScavenger(models, names);
		primary.physicalEdgeMinLen=0;
		check(primary.scavenge("physical", fragment, 0, new ArrayList<int[]>()).isEmpty(),
			"Strict full-model verification must reject the fragment before the physical-edge retry");
		final ArrayList<Orf> ordinary=edgeScavenger(models, names).scavenge("physical", fragment, 0, new ArrayList<int[]>());
		check(ordinary.size()==1 && ordinary.get(0).ncrnaPartial==3 && names[0].equals(ordinary.get(0).trnaModel),
			"The unchanged physical-edge verifier must retain both missing ends and the actual source model");
		final String baseline=ordinary.get(0).toGff();
		for(final int cap:new int[]{50, 384}){
			final NcrnaScavenger caller=edgeScavenger(models, names);
			caller.useDdlProfiles(resource, cap);
			final ArrayList<Orf> calls=caller.scavenge("physical", fragment, 0, new ArrayList<int[]>());
			check(calls.size()==1 && calls.get(0) instanceof NcrnaDdlOrf && calls.get(0).ncrnaPartial==3,
				"Both hybrid caps must preserve the physical partial and its experimental proposal metadata");
			final String observed=calls.get(0).toGff();
			final String prefix=baseline+";ddl_stream=native;ddl_native_shared=";
			check(observed.startsWith(prefix) && Integer.parseInt(observed.substring(prefix.length()))>0,
				"Clipped verification must preserve ordinary GFF bytes and bind positive native support to the winning model");
		}
	}

	private static NcrnaScavenger edgeScavenger(final byte[][] models, final String[] names){
		return edgeScavenger(models, names, true);
	}

	private static NcrnaScavenger edgeScavenger(final byte[][] models, final String[] names, final boolean nativeRanking){
		assert(models.length==2 && names.length==2) : "Edge fixtures bind two independently voted original models";
		final LongHashSet seeds=new LongHashSet();
		final long mask=(1L<<34)-1;
		for(final byte[] model:models){
			long key=0;
			for(int i=0; i<model.length; i++){
				key=((key<<2)|AminoAcid.baseToNumber[model[i]])&mask;
				if(i>=16){seeds.add(key);}
			}
		}
		final NcrnaScavenger caller=new NcrnaScavenger(models, null, names, seeds,
			17, 600, 1500, 9, 384, true, 2f, 0f, 0f, 0, 0f, 1f, 1f, 1f,
			null, null, null, null, -1, -1, -1, -1, 0f,
			NcrnaFamily.LEGACY_START_OFFSETS, NcrnaFamily.LEGACY_STOP_OFFSETS, 0f, 0f, nativeRanking);
		caller.family="cobalamin_rf00174";
		if(nativeRanking){caller.setShortlistPolicy(384, true, 2f, 0f, 0f, 0, 8);}
		else{caller.setModelPatience(8);}
		caller.physicalEdgeMinLen=600;
		caller.scavengePass2=false;
		return caller;
	}

	/** Run last: native blacklist loading has no public undo operation. The suite
	 * exits its own JVM after proving raw scoring neither consumes nor changes it. */
	private static void globalIndependence(final String root, final NcrnaDdlProfileResource resource){
		check(!DynamicDemiLog.blacklistExists(), "This last fixture must establish its own unrelated blacklist");
		final NcrnaDdlProfileResource.Scorer scorer=resource.newScorer();
		final byte[] query=bytes("AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCC");
		final float[] expected=new float[2], observed=new float[2];
		scorer.score(query, expected);
		check(expected[0]>0 && expected[1]>0, "Both literal model profiles must match before changing global DDL state");
		final String path=root+"/unrelated_blacklist.fa";
		try(Output out=new Output(path)){
			out.writer.println(">unrelated_k3\nACGTACGTACGT");
		}
		DynamicDemiLog.setExponent(6);
		DynamicDemiLog.loadBlacklist(path, 3);
		for(final byte[] bases:new byte[][]{query, bytes(""), query}){
			scorer.score(bases, observed);
			check(Arrays.equals(observed, bases.length==0 ? new float[2] : expected),
				"Raw resource scores must replace scratch and remain identical with global exponent6 and an unrelated blacklist");
		}
		final char[] counts={0, 12};
		final NcrnaDdlInterleave first=new NcrnaDdlInterleave(2), second=new NcrnaDdlInterleave(2);
		first.build(counts, expected);
		second.build(counts, observed);
		for(int model=0; model<2; model++){
			check(first.rank(model)==second.rank(model), "Unrelated global state cannot change the raw profile's hybrid order");
		}
		check(DynamicDemiLog.exponentBits()==6 && DynamicDemiLog.blacklistExists(),
			"Caller raw scoring must preserve the process exponent and existing blacklist");
		final DynamicDemiLog nativeSketch=DynamicDemiLog.create(2048, 10, 12345, 0, false, false);
		try{
			new NcrnaDdlQueryKeys(2048, 10).build(query, nativeSketch, false, true);
		}catch(IllegalStateException failure){
			check(failure.getMessage().contains("native sketch destination"), "The construction-only guard must explain its native destination");
			return;
		}
		throw new AssertionError("A supplied native destination must retain its exponent4/no-blacklist construction guard");
	}

	/** Checks the builder CLI on the existing accepted two-member construction
	 * join. These literal names/sequences/votes are independent of its new output. */
	private static void verifyBuilder(final String path){
		final String[] names={"cobalamin_RF00174_fixture0 n=1", "cobalamin_RF00174_fixture1 n=1"};
		final byte[][] sequences={bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("GGGGGGGGGGGGGGGGGGGG")};
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, sequences);
		final NcrnaDdlProfile[] expected=resource.compact() ? compactProfiles(sequences, resource.buckets())
			: resource.dual() ? dualProfiles(sequences) : profiles(sequences, 1);
		for(int model=0; model<2; model++){
			check(resource.members(model)==1, "The CLI must consume each original kept member exactly once");
			for(int bucket=0; bucket<resource.buckets(); bucket++){
				if(resource.key(model, bucket)!=expected[model].keys()[bucket]){
					throw new AssertionError("CLI profile differs from its independent literal member: model="+model+", bucket="+bucket);
				}
			}
		}
		if(resource.compact()){
			compactParity(resource, expected);
		}else if(resource.dual()){
			final float[] scores=new float[2];
			resource.newScorer().score(sequences[0], scores);
			check(Arrays.equals(scores, new float[]{1.5f,0f}), "The dual CLI must preserve both A20 half matches and no G20 match");
		}else{parity(resource, expected);}
	}

	private static NcrnaDdlProfile[] profiles(final byte[][] sequences, final int firstMembers){
		assert(sequences.length==2 && firstMembers>0) : "The literal round-trip fixture has two independently voted positive-member models";
		final NcrnaDdlProfile[] profiles=new NcrnaDdlProfile[2];
		final long wordMask=AbstractKmerTable.makeMiddleMask(10, false, 1);
		for(int model=0; model<profiles.length; model++){
			final DynamicDemiLog sketch=DynamicDemiLog.create(2048, 10, 12345, 0, false, false);
			CardinalityTracker.forEachKmer(sequences[model], 10, 10, false, false,
				(word, length)->sketch.add(word&wordMask));
			final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(2048);
			for(int member=0; member<(model==0 ? firstMembers : 1); member++){votes.add(sketch.maxArray());}
			profiles[model]=votes.finish();
		}
		return profiles;
	}

	/** Bounded literal collision search: two distinct raw keys in one bucket,
	 * repeated query words, shared postings, empty models/query and signed keys. */
	private static void inverted(){
		final NcrnaDdlQueryKeys query=new NcrnaDdlQueryKeys(32, 4);
		final int[] firstWords=new int[32], firstKeys=new int[32];
		Arrays.fill(firstWords, -1);
		int first=-1, second=-1, bucket=-1, key1=-1, key2=-1;
		for(int value=0; value<256 && second<0; value++){
			query.build(word(value), null, false, false);
			check(query.keyCount()==1, "An intact forward4mer must produce one raw query key");
			final int packed=query.packedKey(0), index=packed>>>16, key=packed&65535;
			if(firstWords[index]>=0 && firstKeys[index]!=key){
				first=firstWords[index];
				second=value;
				bucket=index;
				key1=firstKeys[index];
				key2=key;
			}else{
				firstWords[index]=value;
				firstKeys[index]=key;
			}
		}
		check(second>=0, "The fixed256-word fixture must supply distinct raw keys sharing one of32buckets");
		final byte[] sequence=new byte[20];
		Arrays.fill(sequence, (byte)'N');
		System.arraycopy(word(first), 0, sequence, 0, 4);
		System.arraycopy(word(second), 0, sequence, 8, 4);
		System.arraycopy(word(first), 0, sequence, 16, 4);
		query.build(sequence, null, false, false);
		// Native ambiguity-state hashing can retain intermediate keys around N
		// separators. The vector must match its actual unique set,not a two-key
		// assumption about independently hashed blocks.
		final int[] raw=query.packedKeys(), live=new int[query.keyCount()];
		final int[] borrowed=query.packedArray();
		check(raw!=borrowed, "The copying diagnostic accessor must remain distinct from borrowed worker storage");
		for(int i=0; i<live.length; i++){
			live[i]=borrowed[i];
			check(live[i]==query.packedKey(i), "Borrowed storage retains the indexed raw live prefix");
		}
		Arrays.sort(raw);
		Arrays.sort(live);
		System.err.println("INVERTED_FIXTURE_SEQUENCE\t"+new String(sequence, StandardCharsets.US_ASCII));
		System.err.println("INVERTED_FIXTURE_RAW_KEYS\t"+Arrays.toString(raw));
		System.err.println("INVERTED_FIXTURE_LIVE_KEYS\t"+Arrays.toString(live));
		check(Arrays.equals(raw, live), "Inverted iteration must equal the old unique raw-key set,including native ambiguity-state keys; raw="
			+Arrays.toString(raw)+", live="+Arrays.toString(live));
		final char[][] profiles=new char[4][32];
		profiles[0][bucket]=(char)key1;
		profiles[1][bucket]=(char)key2;
		profiles[2][bucket]=(char)key1;
		int absent=65535;
		while(Arrays.binarySearch(raw, (bucket<<16)|absent)>=0){absent--;}
		check(absent>0, "The tiny query cannot occupy all65535nonempty keys in its fixture bucket");
		profiles[3][bucket]=(char)absent;
		final DDLIndexSmall index=new DDLIndexSmall(profiles);
		final int[] scores={-1,-1,-1,-1};
		query.count(index, scores);
		check(Arrays.equals(scores, new int[]{1,1,1,0}), "Distinct keys in one bucket must score only their own shared postings,once per model");
		profiles[0][bucket]=0;
		query.count(index, scores);
		check(scores[0]==1, "An immutable inverted index cannot alias later profile-array mutation");
		profiles[0][bucket]=(char)key1;
		query.count(new DDLIndexSmall(new char[4][32]), scores);
		check(Arrays.equals(scores, new int[4]), "Empty profile buckets must clear all previous worker scores");
		boolean rejected=false;
		try{query.count(new DDLIndexSmall(new char[4][33]), scores);}
		catch(IllegalArgumentException expectedFailure){rejected=true;}
		check(rejected, "The raw collector must reject a shared index with a different width");
		highKeyMatch();
		final int[] reference=new int[4];
		final float[] referenceScores=new float[4], observedScores=new float[4];
		final NcrnaDdlInterleave expected=new NcrnaDdlInterleave(4), observed=new NcrnaDdlInterleave(4);
		final char[] counts={2,0,4,2};
		final java.util.Random random=new java.util.Random(20261113);
		for(int sample=0; sample<64; sample++){
			final byte[] bases=new byte[sample*7];
			for(int i=0; i<bases.length; i++){bases[i]=(byte)"ACGTN".charAt(random.nextInt(5));}
			query.build(bases, null, false, false);
			query.count(index, scores);
			for(int model=0; model<4; model++){
				reference[model]=query.equalBuckets(profiles[model]);
				referenceScores[model]=(float)reference[model];
				observedScores[model]=(float)scores[model];
			}
			check(Arrays.equals(scores, reference), "Inverted scores must exactly equal the old per-model scorer,including empty queries");
			final int before=expected.buildCaller(counts, referenceScores, 2, 4);
			final int after=observed.buildCaller(counts, observedScores, 2, 4);
			check(before==after, "Equal raw scores must retain the native/sketch proposal union size");
			for(int model=0; model<4; model++){
				check(expected.rank(model)==observed.rank(model), "Inversion must preserve every deterministic hybrid model rank");
			}
		}
	}

	/** Exercise a negative packed SEARCH TARGET, not merely an absent high key. */
	private static void highKeyMatch(){
		final NcrnaDdlQueryKeys query=new NcrnaDdlQueryKeys(32, 8);
		int packed=-1;
		for(int value=0; value<65536 && packed<0; value++){
			query.build(word(value, 8), null, false, false);
			check(query.keyCount()==1, "Each complete8mer supplies one real encoded search target");
			final int observed=query.packedKey(0);
			if((observed&32768)!=0){packed=observed;}
		}
		check(packed>=0, "The fixed8mer space must contain a real exponent4key with bit15set to test signed search targets");
		final int bucket=packed>>>16, key=packed&65535;
		final char[][] profiles=new char[3][32];
		profiles[0][bucket]=(char)key;
		profiles[1][bucket]=(char)Math.max(1, key&32767);
		profiles[2][bucket]=(char)(key==65535 ? key-1 : key+1);
		final int[] scores={-1,-1,-1};
		query.count(new DDLIndexSmall(profiles), scores);
		check(Arrays.equals(scores, new int[]{1,0,0}), "A negative packed target must match its own model,not adjacent negative or positive-key postings");
		final float[] hybridScores=new float[scores.length];
		for(int model=0; model<profiles.length; model++){
			check(scores[model]==query.equalBuckets(profiles[model]), "Signed-target lookup must match the old scalar bucket oracle");
			hybridScores[model]=(float)scores[model];
		}
		final NcrnaDdlInterleave order=new NcrnaDdlInterleave(3);
		check(order.buildCaller(new char[3], hybridScores, 2, 3)==1 && order.rank(0)==1
			&& order.rank(1)==0 && order.rank(2)==0,
			"The literal high-key match must be the sole sketch proposal when native support is absent");
	}

	private static byte[] word(final int value){return word(value, 4);}
	private static byte[] word(final int value, final int k){
		check((k==4 || k==8) && value>=0 && value<(1<<(2*k)), "Literal word fixtures enumerate only the declared4mer/8mer spaces");
		final byte[] bases=new byte[k];
		for(int i=0; i<bases.length; i++){bases[i]=(byte)"ACGT".charAt((value>>>(2*(k-1-i)))&3);}
		return bases;
	}

	private static void parity(final NcrnaDdlProfileResource resource, final NcrnaDdlProfile[] profiles){
		final NcrnaDdlProfileResource.Scorer first=resource.newScorer(), second=resource.newScorer();
		final NcrnaDdlRanks original=new NcrnaDdlRanks(profiles, 10, false, true);
		final NcrnaDdlInterleave expected=new NcrnaDdlInterleave(2), observed=new NcrnaDdlInterleave(2);
		final float[] scores=new float[2], other=new float[2], reference=new float[2];
		final char[] counts={0, 12};
		final String[] queries={"AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCC", "aaaaaaaaaaaaaaaaaaaa", "",
			"NNNNNNNNNNNN", "CCCCCCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGGGGGG", "ACGTACGTACGTACGTACGTACGT", "A"};
		for(final String query:queries){
			final byte[] sequence=bytes(query);
			original.score(sequence);
			for(int model=0; model<2; model++){reference[model]=original.score(NcrnaDdlRanks.ALL_KEYS, model);}
			first.score(sequence, scores);
			second.score(bytes("CCCCCCCCCCCCCCCCCCCC"), other);
			first.score(sequence, scores);
			check(Arrays.equals(scores, reference), "Round-trip raw scores must match the completed grid's actual scoring definition");
			expected.build(counts, reference);
			observed.build(counts, scores);
			for(int model=0; model<2; model++){
				check(expected.rank(model)==observed.rank(model), "Resource round-trip must preserve every duplicate-skipping model position");
			}
		}
		check(resource.models()==2, "Independent worker scoring cannot mutate the shared model roster");
	}

	private static void literalKeys(final String root, final String[] names, final byte[][] sequences){
		final char[] keys=new char[2048];
		keys[7]=Character.MAX_VALUE;
		keys[2047]=1;
		final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(2048);
		votes.add(keys);
		final NcrnaDdlProfile profile=votes.finish();
		final String path=root+"/literal_keys.tsv";
		NcrnaDdlProfileResource.write(path, new NcrnaDdlProfile[]{profile, profile}, names, sequences);
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(path, names, sequences);
		check(resource.key(0, 0)==0 && resource.key(0, 7)==65535 && resource.key(1, 2047)==1,
			"Empty, maximum16bit and last-bucket keys round-trip without signed conversion");
	}

	private static void fault(final String path, final String fault, final NcrnaDdlProfile[] profiles,
			final String[] names, final byte[][] sequences){
		assert(profiles.length==2 && names.length==2 && sequences.length==2) : "Each binding fault changes only the literal two-model resource";
		try(Output out=new Output(path)){
			out.writer.println(fault.equals("header") ? NcrnaDdlProfileResource.FORMAT.replace("RAW_ALL_QUERY_KEYS", "SSU") : NcrnaDdlProfileResource.FORMAT);
			final int rows=fault.equals("missing") ? 1 : fault.equals("extra") ? 3 : 2;
			final ByteBuilder b=new ByteBuilder();
			for(int model=0; model<rows; model++){
				final int source=model&1;
				b.clear().append(fault.equals("order") && model==1 ? 0 : model).tab();
				final String name=fault.equals("duplicate_id") ? names[0] : fault.equals("name") && model==0 ? "cobalamin_RF00174_wrong" : names[source];
				b.append(prok.RrnaResourceIO.first(name).getBytes(StandardCharsets.UTF_8)).tab();
				for(int base=0; base<sequences[source].length; base++){
					final byte value=sequences[source][base];
					b.append(fault.equals("sequence") && model==0 && base==0 ? (byte)'G' : value>='a' && value<='z' ? (byte)(value-32) : value);
				}
				b.tab().append(fault.equals("members") && model==0 ? 0 : profiles[source].members);
				final int buckets=fault.equals("short") && model==0 ? 2047 : 2048;
				for(int bucket=0; bucket<buckets; bucket++){
					final int value=model==0 && bucket==7 && fault.equals("negative") ? -1
						: model==0 && bucket==7 && fault.equals("overflow") ? 65536 : profiles[source].keys()[bucket];
					b.tab().append(value);
				}
				out.writer.print(b.nl());
			}
		}
	}
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void reject(final Runnable action, final String reason){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(reason), "Wrong resource rejection: "+failure.getMessage());
			return;
		}
		throw new AssertionError("Expected resource rejection containing: "+reason);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
