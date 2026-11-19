package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GeneCaller;
import prok.GeneModel;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaDdlShortlist;
import prok.NcrnaDdlOrf;
import prok.NcrnaFamily;
import prok.NcrnaAlignmentWork;
import prok.NcrnaCpuProfile;
import prok.NcrnaModelAttemptInstrumentSink;
import prok.NcrnaScavenger;
import prok.NcrnaWindowCaptureSink;
import prok.TrnaConsensusBuilder;
import prok.Orf;
import shared.Shared;
import stream.Read;

/** Finite real-verifier test: the exact model lies beyond the finite384-model cap.
 * Zero profile keys deliberately force ordinal ties, independently of voting.
 * This is caller plumbing and observation, not a biological or timing result.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlUnlimitedTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlUnlimitedTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown caller fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory without input, overwrite, append or sampling");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "The finite caller fixture requires assertions");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture output must be a new directory");
			DynamicDemiLog.setExponent(5);
			run(root.getPath());
			rescue(root.getPath());
			operands(root.getPath());
			spanEligibility(root.getPath());
			check(DynamicDemiLog.exponentBits()==5, "Selected caller queries must preserve CallGenes' normal global exponent5");
			pp.outstream.println("NCRNA_DDL_UNLIMITED_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Real Caller          ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final String root){
		final int count=400, length=320;
		final byte[] query=Prok16sPartialRescueTest.random(length, 20261115);
		final byte[][] sequences=new byte[count][];
		final String[] names=new String[count];
		final int[] voters=new int[count];
		Arrays.fill(voters, 1);
		for(int model=0; model<count; model++){
			sequences[model]=query.clone();
			names[model]="cobalamin_RF00174_unlimited_"+model;
			if(model==0){sequences[model][length/2]=different(query[length/2]);}
			else if(model<count-1){
				for(int base=0; base<length; base+=4){sequences[model][base]=different(query[base]);}
			}
		}
		final String path=root+"/ordinal_ties.tsv";
		NcrnaDdlProfileResource.writeSelected(path, new char[count][512], voters, names, sequences, 512, new LongHashSet(4));
		final NcrnaDdlProfileResource profiles=NcrnaDdlProfileResource.load(path, names, sequences);
		final NcrnaDdlShortlist order=new NcrnaDdlShortlist(profiles, Integer.MAX_VALUE);
		final int[] first=order.build(query, null);
		check(order.unlimited() && order.size()==count && first[count-1]==count-1,
			"An unlimited zero-score order must include every model, including the exact model after position384");
		check(order.build(new byte[0], null)==first && order.size()==count,
			"Empty queries retain the complete zero-score order while reusing caller scratch");
		final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			final LongHashSet seeds=new LongHashSet(4);
			long seed=0;
			for(int base=0; base<17; base++){seed=(seed<<2)|dna.AminoAcid.baseToNumber[query[base]];}
			seeds.add(seed);
			final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", sequences, null, names, seeds,
				17, 109, length, 9, Integer.MAX_VALUE, true, 2f, 0f, 0f, 0, 0f, 20f, .99f, .90f);
			family.modelPatience=Integer.MAX_VALUE;
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(family);
			final GeneCaller caller=new GeneCaller(60, 10, 10, -1f, -1f, -1f, -1f, -1f, new GeneModel(false));
			caller.setCobalaminDdl(profiles, Integer.MAX_VALUE);
			final NcrnaScavenger unlimited=caller.ncrnaScavenger(0);
			unlimited.trimAlignmentExtent=false;
			unlimited.scavengePass2=false;
			check(!unlimited.hasNativeRankingIndex() && unlimited.indexTopN==Integer.MAX_VALUE
				&& unlimited.modelPatience()==Integer.MAX_VALUE,
				"The actual GeneCaller factory must consume both unlimited controls without constructing a native index");
			final int[] attempts={0}, windows={0};
			unlimited.setWindowCaptureSink((contig, strand, pass, start, stop, source, bases, window)->{
				check(start==0 && stop==length-1 && Arrays.equals(window, query),
					"Capture must observe the actual full candidate window before ranking");
				windows[0]++;
			});
			unlimited.setModelAttemptSink((contig, strand, pass, start, stop, rank, model, shared, identity, span, passed)->{
				check(shared==-1 && rank==model+1 && rank==attempts[0]+1,
					"Each actual primary attempt must retain its ordinal and mark unmeasured native support as unavailable");
				check(unlimited.ddlScore(model)==0f, "Every observed score must come from this zero-key current-window order");
				if(model==count-1){check(identity==1f && passed, "The final exact model must genuinely pass the real aligner");}
				attempts[0]++;
			});
			final ArrayList<Orf> calls=measured(unlimited, query, "unlimited", count);
			check(calls.size()==1 && names[count-1].equals(calls.get(0).trnaModel)
				&& attempts[0]==count && unlimited.alignmentCount()==count && windows[0]==1,
				"All400 real verifications must select the exact last model despite an earlier passing model, with one captured window");
			final NcrnaDdlOrf best=(NcrnaDdlOrf)calls.get(0);
			check(best.acceptedIdentity()==1f && best.proposalRank()==count && best.proposalScore()==0f,
				"The accepted winner must retain its actual identity, zero DDL score and late order without inference from path score");
			check(best.toGff().contains(";ddl_stream=sketch_only;ddl_native_shared=NA"),
				"Final call provenance must retain unavailable native support after winning-model selection");
			final NcrnaFamily finiteFamily=new NcrnaFamily("cobalamin_rf00174", sequences, null, names, seeds,
				17, 109, length, 9, 384, true, 2f, 0f, 0f, 0, 0f, 20f, .99f, .90f);
			finiteFamily.modelPatience=8;
			GeneCaller.ncrnaFamilies.set(0, finiteFamily);
			final GeneCaller finiteCaller=new GeneCaller(60, 10, 10, -1f, -1f, -1f, -1f, -1f, new GeneModel(false));
			finiteCaller.setCobalaminDdl(profiles, 384);
			final NcrnaScavenger finite=finiteCaller.ncrnaScavenger(0);
			finite.trimAlignmentExtent=false;
			finite.scavengePass2=false;
			final ArrayList<Orf> earlier=measured(finite, query, "finite", 9);
			check(earlier.size()==1 && names[0].equals(earlier.get(0).trnaModel) && finite.alignmentCount()<384,
				"The explicit finite control must retain its early passing model and patience stopping, proving the all-model effect");
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.addAll(saved);
		}
	}

	/** Both candidates pass HBM, while the earlier one has the higher HBM score
	 * and the later one has the higher identity. This distinguishes finite legacy
	 * HBM selection from the explicit unlimited best-accepted-identity baseline. */
	private static void rescue(final String root){
		final byte[] query=Prok16sPartialRescueTest.random(96, 71115);
		final byte[][] sequences={query.clone(), query.clone()};
		sequences[0][32]=different(query[32]);
		sequences[0][64]=different(query[64]);
		sequences[1][32]=different(query[32]);
		final String[] names={"cobalamin_RF00174_rescue0", "cobalamin_RF00174_rescue1"};
		final BaseGraph[] graphs=new BaseGraph[2];
		final float[] identities=new float[2], hbms=new float[2];
		for(int model=0; model<2; model++){
			graphs[model]=new BaseGraph(names[model], sequences[model], null, model, 0);
			if(model==0){
				final Read support=new Read(query.clone(), null, "support", 0);
				support.start=0;
				support.stop=query.length-1;
				support.match=new byte[query.length];
				for(int base=0; base<query.length; base++){support.match[base]=query[base]==sequences[model][base] ? (byte)'m' : (byte)'S';}
				for(int copy=0; copy<100; copy++){graphs[model].add(support);}
			}
			BaseGraphHelper.initForScoring(graphs[model]);
			final AlignmentStats trace=new AlignmentStats(true);
			identities[model]=QuantumAligner.alignAndTraceStatic(sequences[model], query, trace);
			modelRead(query, sequences[model], trace);
			hbms[model]=TrnaConsensusBuilder.scoreAlignedAgainstModel(query, trace, graphs[model]);
		}
		check(identities[1]>identities[0] && identities[1]<1f && hbms[0]>hbms[1] && hbms[1]>=.5f,
			"The real rescue fixture must distinguish identity from HBM rank, with both models below direct acceptance and above HBM acceptance");
		final String path=root+"/rescue.tsv";
		NcrnaDdlProfileResource.writeSelected(path, new char[2][512], new int[]{1,1}, names, sequences, 512, new LongHashSet(4));
		final NcrnaDdlProfileResource profiles=NcrnaDdlProfileResource.load(path, names, sequences);
		final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			for(final int cap:new int[]{384,Integer.MAX_VALUE}){
				final NcrnaFamily family=family(sequences, graphs, names, query, cap, 64, 1f, .9f, .5f);
				final NcrnaScavenger worker=worker(family, profiles, cap);
				final int expected=cap==384 ? 0 : 1;
				final int[] phases={0,0}, completed={0};
				worker.setModelAttemptSink(new NcrnaModelAttemptInstrumentSink(){
					@Override
					public void modelAttempt(String name, int strand, int pass, int start, int stop,
							int rank, int model, int shared, float identity, int span, boolean accepted){
						check(shared==-1 && !accepted && rank==model+1 && identity==identities[model],
							"Short-window primary attempts must use the declared Quantum operands and remain below direct acceptance");
						phases[0]++;
					}
					@Override
					public void rescueAttempt(String name, int strand, int pass, int start, int stop,
							int candidateStart, int candidateStop, int rank, int model, float identity,
							int span, float hbm, boolean hbmMeasured, boolean accepted){
						check(candidateStart==0 && candidateStop==query.length-1 && span==query.length
							&& rank==model+1 && identity==identities[model] && hbmMeasured && hbm==hbms[model] && accepted,
							"Each rescue must reuse its one Quantum trace for the actual reference span, identity and HBM score");
						phases[1]++;
					}
				});
				worker.setWindowCaptureSink(new NcrnaWindowCaptureSink(){
					@Override
					public void capture(String name, int strand, int pass, int start, int stop,
							int source, byte[] bases, byte[] window){
						check(Arrays.equals(window, query), "Rescue must retain the actual primary window");
					}
					@Override
					public void completed(String name, int strand, int pass, int start, int stop, Orf accepted, long alignments){
						check(alignments==4 && accepted instanceof NcrnaDdlOrf
							&& ((NcrnaDdlOrf)accepted).acceptedIdentity()==identities[expected],
							"Completed-window accounting must include both primary and both rescue kernels and retain the actual winner");
						completed[0]++;
					}
				});
				final ArrayList<Orf> calls=measured(worker, query, "rescue", 4, 2);
				check(calls.size()==1 && names[expected].equals(calls.get(0).trnaModel)
					&& phases[0]==2 && phases[1]==2 && completed[0]==1 && worker.hbmScoreCalls()==2,
					"Unlimited must choose later higher accepted identity; finite retains higher HBM, with all four kernels and two HBM evaluations conserved");
			}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.addAll(saved);
		}
	}

	/** Model query is shorter, equal and longer than its genomic reference.
	 * The third case has a real internal deletion, not an operand swap. */
	private static void operands(final String root){
		final byte[] core=Prok16sPartialRescueTest.random(96, 8193);
		for(int mode=0; mode<3; mode++){
			final byte[] query=core.clone();
			final byte[] model;
			final int expectedStart=mode==0 ? 16 : 0, expectedStop=mode==0 ? 79 : 95;
			if(mode==0){model=Arrays.copyOfRange(core, 16, 80);}
			else if(mode==1){model=core.clone();}
			else{
				model=new byte[112];
				System.arraycopy(core, 0, model, 0, 48);
				System.arraycopy(Prok16sPartialRescueTest.random(16, 96112), 0, model, 48, 16);
				System.arraycopy(core, 48, model, 64, 48);
			}
			final AlignmentStats trace=new AlignmentStats(true);
			final float expectedIdentity=QuantumAligner.alignAndTraceStatic(model, query, trace);
			final float primaryIdentity=QuantumAligner.alignStatic(model, query, new int[4]);
			final float directPass=.75f;
			check(primaryIdentity>=directPass, "The operand fixture must exercise direct acceptance before optional trimming: mode="
				+mode+" primary="+primaryIdentity+" threshold="+directPass);
			check(trace.qLen==model.length && trace.rLen==query.length
				&& trace.rStart==expectedStart && trace.rStop==expectedStop,
				"Quantum positions must remain genomic-reference coordinates for every model/window length relation");
			modelRead(query, model, trace);
			final String[] names={"cobalamin_RF00174_operand_"+mode};
			final byte[][] sequences={model};
			final String path=root+"/operand_"+mode+".tsv";
			NcrnaDdlProfileResource.writeSelected(path, new char[1][512], new int[]{1}, names, sequences, 512, new LongHashSet(4));
			final NcrnaDdlProfileResource profiles=NcrnaDdlProfileResource.load(path, names, sequences);
			final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
			try{
				for(final boolean trim:new boolean[]{false,true}){
					final NcrnaFamily family=family(sequences, null, names, query, Integer.MAX_VALUE, 32, directPass, .7f);
					family.trimAlignmentExtent=trim;
					final NcrnaScavenger worker=worker(family, profiles, Integer.MAX_VALUE);
					final String scenario="operands_mode_"+mode+"_trim_"+trim+"_primary_"+primaryIdentity
						+"_direct_"+(primaryIdentity>=directPass)+"_trace_"+expectedIdentity;
					final ArrayList<Orf> calls=measured(worker, query, scenario, trim ? 2 : 1);
					check(calls.size()==1 && calls.get(0).start==expectedStart && calls.get(0).stop==expectedStop
						&& ((NcrnaDdlOrf)calls.get(0)).acceptedIdentity()==expectedIdentity,
						"Selected short-window calling and optional Quantum trimming must preserve the actual genomic span and accepted identity");
				}
			}finally{
				GeneCaller.ncrnaFamilies.clear();
				GeneCaller.ncrnaFamilies.addAll(saved);
			}
		}
	}

	/** A perfect short decoy must not hide a longer passing gene in either
	 * selected walk; the global gene floor is part of eligibility, not ranking. */
	private static void spanEligibility(final String root){
		final byte[] query=Prok16sPartialRescueTest.random(128, 96128);
		final byte[][] sequences={Arrays.copyOfRange(query, 32, 96), query.clone()};
		for(int base=3; base<query.length; base+=8){sequences[1][base]=different(query[base]);}
		final String[] names={"cobalamin_RF00174_short_decoy", "cobalamin_RF00174_long_pass"};
		final String path=root+"/span_eligibility.tsv";
		NcrnaDdlProfileResource.writeSelected(path, new char[2][512], new int[]{1,1}, names, sequences, 512, new LongHashSet(4));
		final NcrnaDdlProfileResource profiles=NcrnaDdlProfileResource.load(path, names, sequences);
		final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			for(final int cap:new int[]{384,Integer.MAX_VALUE}){
				for(final float pass:new float[]{.8f,.95f}){
					final NcrnaFamily family=family(sequences, null, names, query, cap, 96, pass, .7f);
					final NcrnaScavenger worker=worker(family, profiles, cap);
					final ArrayList<Orf> calls=measured(worker, query, "span_eligibility", pass==.8f ? 2 : 4);
					if(pass==.8f){
						check(calls.size()==1 && names[1].equals(calls.get(0).trnaModel)
							&& calls.get(0).stop-calls.get(0).start+1>=96,
							"A perfect64-base model below minLen96 cannot hide the later128-base passing model in finite or unlimited selected verification");
					}else{
						check(calls.isEmpty(), "When the longer model fails identity, rescue cannot use a perfect sub-minimum fragment to certify the longer candidate interval");
					}
				}
			}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.addAll(saved);
		}
	}

	private static Read modelRead(final byte[] query, final byte[] model, final AlignmentStats trace){
		final Read aligned=TrnaConsensusBuilder.invertToModelFrame(query, trace);
		check(aligned!=null && aligned.start==0 && aligned.stop==model.length-1
			&& consumption(aligned.match, true)==aligned.length()
			&& consumption(aligned.match, false)==model.length,
			"Consensus-as-query traceback must invert to an HBM Read consuming the actual model and candidate span, without swapping coordinate roles");
		return aligned;
	}

	private static ArrayList<Orf> measured(final NcrnaScavenger worker, final byte[] query, final String name, final long expected){
		return measured(worker, query, name, expected, 0);
	}

	private static ArrayList<Orf> measured(final NcrnaScavenger worker, final byte[] query, final String name,
			final long expected, final long expectedHbm){
		final NcrnaCpuProfile profile=NcrnaCpuProfile.create(1);
		check(profile!=null, "Use the fixture launcher so native kernel accounting is enabled before class initialization");
		final long before=worker.alignmentCount();
		final ArrayList<Orf> calls;
		profile.begin(0, NcrnaCpuProfile.SCAVENGE);
		try{calls=worker.scavenge(name, query, 0, new ArrayList<int[]>());}
		finally{profile.end();}
		final NcrnaAlignmentWork work=profile.alignmentWork(0);
		check(work.logical==expected && work.physical==expected && work.hits==0 && work.hbmCalls==expectedHbm
			&& work.windowPhysical==expected && work.windowLogical==expected && work.windows==1
			&& worker.alignmentCount()-before==expected,
			"Native physical kernels, logical requests, whole-window totals and caller counts must reconcile, including rescue and endpoint work: "+name
			+" expected="+expected+" logical="+work.logical+" physical="+work.physical+" hits="+work.hits
			+" expectedHBM="+expectedHbm+" hbm="+work.hbmCalls+" windowPhysical="+work.windowPhysical
			+" windowLogical="+work.windowLogical+" windows="+work.windows+" caller="+(worker.alignmentCount()-before));
		return calls;
	}

	private static NcrnaFamily family(final byte[][] sequences, final BaseGraph[] graphs, final String[] names,
			final byte[] query, final int cap, final int minimum, final float pass, final float borderline){
		return family(sequences, graphs, names, query, cap, minimum, pass, borderline, .75f);
	}

	private static NcrnaFamily family(final byte[][] sequences, final BaseGraph[] graphs, final String[] names,
			final byte[] query, final int cap, final int minimum, final float pass, final float borderline, final float hbm){
		final LongHashSet seeds=new LongHashSet(4);
		long seed=0;
		for(int base=0; base<17; base++){seed=(seed<<2)|dna.AminoAcid.baseToNumber[query[base]];}
		seeds.add(seed);
		final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", sequences, graphs, names, seeds,
			17, minimum, query.length*2, 9, cap, true, 2f, 0f, 0f, 0, 0f, 20f, pass, borderline,
			hbm, .9f, null, null, null, null, -1, -1, -1, -1, 0f);
		family.modelPatience=cap==Integer.MAX_VALUE ? Integer.MAX_VALUE : 8;
		family.trimAlignmentExtent=false;
		family.scavengePass2=false;
		return family;
	}

	private static NcrnaScavenger worker(final NcrnaFamily family, final NcrnaDdlProfileResource profiles, final int cap){
		GeneCaller.ncrnaFamilies.clear();
		GeneCaller.ncrnaFamilies.add(family);
		final GeneCaller caller=new GeneCaller(60, 10, 10, -1f, -1f, -1f, -1f, -1f, new GeneModel(false));
		caller.setCobalaminDdl(profiles, cap);
		return caller.ncrnaScavenger(0);
	}

	private static int consumption(final byte[] match, final boolean query){
		check(match!=null, "HBM coordinate validation requires the actual normalized traceback");
		int count=0;
		for(final byte op:match){
			check(op=='m' || op=='S' || op=='N' || op=='I' || op=='D', "Unexpected HBM traceback operation");
			if(op!=(query ? 'D' : 'I')){count++;}
		}
		return count;
	}

	private static byte different(final byte base){return base=='A' ? (byte)'C' : (byte)'A';}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
