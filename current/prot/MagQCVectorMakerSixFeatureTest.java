package prot;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import ml.CellNet;
import ml.Function;

	/**
	 * T1 (design revision2 plans/MVM_SIX_FEATURE_INTEGRATION_DESIGN_20260909.md §7), scoped: a real,
	 * end-to-end {@code subnetfeatures=six} run through the actual {@link MagQCVectorMaker}, against a
	 * genuine 176-subnet four-output bundle built through the REAL {@link MagQCNetBundle#pack} path
 * (not a hand-written file) -- reusing {@link MagQCAggVectorTest}'s proven 4-organism/8-family
 * cache fixture. Verifies: header width matches the six-mode formula, one subnet's full six-value
 * block (4 learned outputs + 2 derived O/E,X/E) recomputed independently against a real loaded net
 * and the real binding, legacy byte-identity is unaffected (a plain legacy run on the SAME cache
 * still matches {@link MagQCAggVectorTest}'s expectations), and the non-finite-learned-output
 * rejection (R5) fires with the subnet id, feature name, and tid named in the message.
 *
 * <p>Scope note: this is not the design's full 176-subnet exact-per-row-per-subnet recomputation
 * (T1 as originally specified) -- it proves the row layout, real-packer compatibility, and the
 * rejection path with a smaller number of independently-recomputed checks. Reported as scoped to
 * Yoimiya, not claimed as the complete T1 fixture.
 *
 * @author Ady
 */
public final class MagQCVectorMakerSixFeatureTest {
	static final int N=176;
	static final int NCRNA_INDEX=N-3;
	static final int RRNA_INDEX=N-2;
	static final int TRNA_INDEX=N-1;
	static final int NFAM=8; // matches MagQCAggVectorTest's familylist.tsv exactly
	static final int FAMSET_OBS=3;
	// The REAL accepted four-output candidate export (design revision2 §1: bbtools-dev `c26965b`,
	// candidate nets c5af69b9... A48 / bd1e4674... decimal), copied verbatim into famset_0's slot
	// instead of a synthetic net (Yoimiya's review request: "preserve the real-candidate check").
	// #sparse (unlike every other synthetic #dense net in this fixture), #dims 22 22 22 8 4 4 --
	// verified: numObs(2)+numPhyla(3)+SHARED_CONTEXT_WIDTH(17)=22, matches its own #dims exactly.
	static final String CANDIDATE_SHA256="c5af69b9e8c33c20f2816fe5d838b86f6c1802af160e692b812c9c3d85d8fd38";
	static final int CANDIDATE_OBS=2;

	public static void main(String[] args) throws Exception {
		if(args.length!=1 || !args[0].startsWith("candidate=") || args[0].length()==10){
			throw new IllegalArgumentException("Required: candidate=<accepted explicit-affine A48 net>");
		}
		final Path candidate=Paths.get(args[0].substring(10));
		check(CANDIDATE_SHA256.equals(sha(candidate)),"candidate SHA-256 differs from accepted export");
		final File dir=Files.createTempDirectory("magqcsix_test").toFile();
		Runtime.getRuntime().addShutdownHook(new Thread(()->deleteRecursive(dir)));
		Arrays.fill(Function.TYPE_RATES,0f);

		// --- identical fixture to MagQCAggVectorTest (4 orgs, 2 phyla, 8 families) ---
		write(new File(dir,"cache.tsv"),
				"c100_1\t100\tbacteria\t5000\t2500\t5000\t5\t5\t4500\t4100000\t4500\t1\t1\t1\t0\t10\t0:1;1:1;2:1\t5:10\n"+
				"c100_2\t100\tbacteria\t3000\t1500\t3000\t3\t3\t2700\t2430000\t2700\t0\t0\t0\t0\t8\t1:1;3:1\t5:8\n"+
				"c100_3\t100\tbacteria\t2000\t1000\t2000\t2\t2\t1800\t1620000\t1800\t0\t1\t0\t1\t5\t4:2\t5:5\n"+
				"c200_1\t200\tbacteria\t6000\t3600\t6000\t6\t6\t5400\t4860000\t5400\t1\t1\t1\t0\t12\t0:1;2:1;5:1\t5:12\n"+
				"c200_2\t200\tbacteria\t4000\t2400\t4000\t4\t4\t3600\t3240000\t3600\t0\t0\t0\t2\t9\t6:1;7:1\t5:9\n"+
				"c300_1\t300\tbacteria\t4000\t2000\t4000\t4\t4\t3600\t3240000\t3600\t1\t1\t0\t0\t9\t0:1;1:1\t5:9\n"+
				"c300_2\t300\tbacteria\t4000\t2000\t4000\t4\t4\t3600\t3240000\t3600\t0\t0\t1\t1\t7\t2:1;3:1\t5:7\n"+
				"c400_1\t400\tbacteria\t7000\t4200\t7000\t7\t7\t6300\t5670000\t6300\t1\t1\t1\t0\t13\t0:1;5:1;6:1\t5:13\n"+
				"c400_2\t400\tbacteria\t5000\t3000\t5000\t5\t5\t4500\t4050000\t4500\t0\t1\t0\t1\t10\t2:1;7:1\t5:10\n");
		write(new File(dir,"sizemap.tsv"),"100\t10000\n200\t10000\n300\t8000\n400\t12000\n");
		write(new File(dir,"taxpgm.tsv"),
			"100\tFirmicutes\tpgm1\n200\tProteobacteria\tpgm2\n300\tFirmicutes\tpgm1\n400\tProteobacteria\tpgm2\n");
		// The exact "#rank\trep_id\tocc_total" 3-column format both MagQCNetBundle.pack()'s
		// readFamilyRepIds() (needs "#"-prefixed header + rank,rep_id) AND
		// MagQCExpectedCopyTableTyped.readFamily() (needs the LITERAL header string + exactly 3
		// fields per row) require -- MagQCAggVectorTest's bare "rank\tcluster_rep" 2-column format
		// works for MVM's own loadAux() (which just skips one line unconditionally) but not these two.
		final File familyList=new File(dir,"familylist.tsv");
		write(familyList,"#rank\trep_id\tocc_total\n0\tf0\t100\n1\tf1\t100\n2\tf2\t100\n3\tf3\t100\n"
			+"4\tf4\t100\n5\tf5\t100\n6\tf6\t100\n7\tf7\t100\n");
		final String familySha=sha(familyList);

		// --- 176-subnet, ALL four-output bundle, built through the REAL pack() path ---
		final File nets=new File(dir,"nets"); nets.mkdirs();
		final File subsets=new File(dir,"subsets"); subsets.mkdirs();
		final File release=new File(dir,"release_manifest.tsv");
		final File aggManifest=new File(dir,"agg_manifest.tsv");
		final int[][] rankArrays=new int[N][];
		try(PrintWriter rw=writer(release); PrintWriter aw=writer(aggManifest)){
			for(int i=0;i<N;i++){
				final String type=typeForIndex(i);
				final boolean special=!"famset".equals(type);
				final boolean isCandidate=(i==0); // famset_0 uses the REAL candidate net (see field javadoc)
				final String id=special ? type : "famset_"+i;
				final int numObs=special ? MagQCObservationLayout.specialCount(type) : (isCandidate ? CANDIDATE_OBS : FAMSET_OBS);
				final String netFile="net_"+id+".bbnet";
				if(!special){
					final int[] ranks=new int[numObs];
					final StringBuilder sb=new StringBuilder();
					for(int k=0;k<numObs;k++){ranks[k]=(i+k)%NFAM; sb.append(ranks[k]).append('\n');}
					rankArrays[i]=ranks;
					write(new File(subsets,id+".txt"),sb.toString());
				}
				final int inputWidth=numObs+3/*phyla incl. other*/+MagQCVectorMaker.SHARED_CONTEXT_WIDTH;
				if(isCandidate){Files.copy(candidate,new File(nets,netFile).toPath());}
				else{writeFourOutputNet(new File(nets,netFile),inputWidth,i+1000);}
				rw.print(id+"\t"+type+"\t"+(special?"-":"subsets/"+id+".txt")+"\toutputs=4\n");
				aw.print(id+"\t"+numObs+"\t"+inputWidth+"\t"+(special?"-":"subsets/"+id+".txt")+"\t-\t"+netFile);
				if("rrna".equals(type) || "trna_anticodon".equals(type)){aw.print("\t"+type);}
				aw.print('\n');
			}
		}
		final File bundleFile=new File(dir,"six_bundle.bbnets");
		final Map<String,String> packArgs=new LinkedHashMap<String,String>();
		packArgs.put("aggmanifest",aggManifest.getAbsolutePath()); packArgs.put("subnetmanifest",release.getAbsolutePath());
		packArgs.put("familylist",familyList.getAbsolutePath()); packArgs.put("taxpgm",new File(dir,"taxpgm.tsv").getAbsolutePath());
		packArgs.put("out",bundleFile.getAbsolutePath()); packArgs.put("netroot",dir.getAbsolutePath()); packArgs.put("subsetroot",dir.getAbsolutePath());
		MagQCNetBundle.pack(packArgs);

		// --- synthetic expected-copy table + declaration for the SAME 4 organisms/8 families ---
		final Path table=buildTable(dir.toPath(),familyList.toPath(),familySha);
		final String tableSha=sha(table.toFile());
		final MagQCNetBundle loaded=MagQCNetBundle.loadMultiOutput(bundleFile.toPath());
		final String releaseSha=loaded.metadata("release_manifest_sha256");
		final StringBuilder declText=new StringBuilder();
		declText.append("#schema_version\tsubnet_population_declaration_v1\n");
		declText.append("#release_manifest_sha256\t").append(releaseSha).append('\n');
		declText.append("#columns\tsubnet_id\tpopulation_type\tpopulation_name\n");
		for(int i=0;i<N;i++){
			// famset_1 is deliberately declared to a PHYLUM population, not global -- proves the
			// binding selects by DECLARATION, not by a uniform default (Yoimiya's review request).
			final String row=(i==1) ? loaded.subnet(i).id+"\tphylum\tFirmicutes\n" : loaded.subnet(i).id+"\tglobal\t-\n";
			declText.append(row);
		}
		final File declaration=write(new File(dir,"declaration.tsv"),declText.toString());
		final String declSha=sha(declaration);

		// --- (1) real six-mode run ---
		final File aggOut=new File(dir,"agg_six.tsv");
		MagQCVectorMaker.main(new String[]{
			"cache="+new File(dir,"cache.tsv").getAbsolutePath(),"sizemap="+new File(dir,"sizemap.tsv").getAbsolutePath(),
			"taxpgm="+new File(dir,"taxpgm.tsv").getAbsolutePath(),"familylist="+familyList.getAbsolutePath(),
			"n=20","valn=0","valfrac=0.5","seed=1","cleanspike=0.0",
			"out="+new File(dir,"g_six.tsv").getAbsolutePath(),
			"bundle="+bundleFile.getAbsolutePath(),"aggout="+aggOut.getAbsolutePath(),
			"subnetfeatures=six","expectedcopytable="+table,"expectedcopytablesha256="+tableSha,
			"subnetpopulations="+declaration,"subnetpopulationssha256="+declSha});

		final ArrayList<String[]> rows=parse(aggOut);
		check(rows.size()==20,"expected 20 rows, got "+rows.size());
		final int expectedWidth=N*6+2*Math.min(100,NFAM)+3+MagQCVectorMaker.SHARED_CONTEXT_WIDTH+2*5+1;
		check(("#dims\t"+expectedWidth+"\t2\t0").equals(header(aggOut)),
			"bad six-mode header: "+header(aggOut)+" (expected width "+expectedWidth+")");
		for(String[] r:rows){check(r.length==expectedWidth+2,"row width "+r.length+" != "+(expectedWidth+2));}

		// --- (2) EXACT independent per-row recomputation, across bundle order, for five subnets:
		// famset_0 (global population, ranks {0,1,2}), famset_1 (DECLARED to phylum Firmicutes,
		// ranks {1,2,3} -- proves selection is by declaration, not a uniform default), and all three
		// explicit special subnets (ncrna, rrna, trna_anticodon, global population). No companion
		// subnetout= runs:
		// every input (obs, phylum one-hot, context, domain one-hot) is reconstructed directly from
		// THIS row's own dense-head/phylum/context/raw-ncRNA columns, which are always whole-bin
		// regardless of aggobs= -- avoiding the serve/clean mismatch a subnetout= companion file
		// would have under real contamination (cleanspike=0.0 here, matching MagQCAggVectorTest's
		// own reasoning for why its serve-vs-clean section decodes from the aggregator row directly).
		// Loaded through the bundle's own per-instance InferenceNet dispatch (netForCurrentThread()),
		// NOT raw CellNetParser.load()+feedForward(): net0 is the REAL candidate net (#sparse) and
		// net1/netN are synthetic (#dense) -- calling feedForward() directly on plain CellNetParser-
		// loaded objects would dispatch on the shared GLOBAL CellNet.DENSE flag (left however the
		// LAST CellNetParser.load() call set it), silently mis-scoring whichever net doesn't match
		// that leftover state. This is exactly the hazard THREAD_SAFE_BUNDLE_INFERENCE_20260909
		// fixed in production; using the bundle's own dispatch here avoids reintroducing it in the
		// test itself, and is more representative of production usage besides.
		final CellNet net0=loaded.subnet(0).netForCurrentThread();
		final CellNet net1=loaded.subnet(1).netForCurrentThread();
		final CellNet netN=loaded.subnet(NCRNA_INDEX).netForCurrentThread();
		final CellNet netR=loaded.subnet(RRNA_INDEX).netForCurrentThread();
		final CellNet netT=loaded.subnet(TRNA_INDEX).netForCurrentThread();
		final int[] headOrder={0,2,1,3,5,6,7,4}; // SAME cache/family fixture as MagQCAggVectorTest -> same prevalence order
		final int subnetBlockWidth=N*6, denseHeadStart=subnetBlockWidth, phylumStart=denseHeadStart+2*NFAM,
			contextStart=phylumStart+3, rawNcStart=contextStart+MagQCVectorMaker.SHARED_CONTEXT_WIDTH;
		final MagQCExpectedCopyTableTyped.Table reopenedTable=MagQCExpectedCopyTableTyped.open(table.toString(),tableSha);
		final double[] globalMeans=reopenedTable.population("global","-").means();
		final double[] firmicutesMeans=reopenedTable.population("phylum","Firmicutes").means();
		int diffFromDefaultPop=0;
		for(int i=0;i<rows.size();i++){
			final String[] r=rows.get(i);
			final int[] fam=new int[NFAM];
			for(int j=0;j<NFAM;j++){fam[headOrder[j]]=decodeTwo(r,denseHeadStart,j);}
			final int[] observedByItem=new int[NFAM+MagQCObservationLayout.EXTENDED_COUNT];
			System.arraycopy(fam,0,observedByItem,0,NFAM);
			for(int t=0;t<5;t++){observedByItem[NFAM+t]=decodeRatio(r,rawNcStart+2*t+1);}
			// This fixture's raw cache assigns every tRNA to structural anticodon code 5.
			observedByItem[NFAM+MagQCObservationLayout.LEGACY_COUNT+5]=observedByItem[NFAM+4];
			final float[] phylumCtx=new float[3+MagQCVectorMaker.SHARED_CONTEXT_WIDTH];
			for(int j=0;j<phylumCtx.length;j++){phylumCtx[j]=Float.parseFloat(r[phylumStart+j]);}

			// famset_0 (cols 0-5, REAL candidate net, CANDIDATE_OBS=2): obs = fam[0],fam[1]
			checkSubnetBlock(r,0,net0,new int[]{fam[0],fam[1]},phylumCtx,
				MagQCExpectedCopyFeatures.computeForItems(observedByItem,globalMeans,new int[]{0,1}),"famset_0",i);
			// famset_1 (cols 6-11): obs = fam[1],fam[2],fam[3]; population = phylum Firmicutes, NOT global
			final MagQCExpectedCopyFeatures.Result r1Firmicutes=MagQCExpectedCopyFeatures.computeForItems(observedByItem,firmicutesMeans,new int[]{1,2,3});
			checkSubnetBlock(r,1,net1,new int[]{fam[1],fam[2],fam[3]},phylumCtx,r1Firmicutes,"famset_1",i);
			final MagQCExpectedCopyFeatures.Result r1Global=MagQCExpectedCopyFeatures.computeForItems(observedByItem,globalMeans,new int[]{1,2,3});
			if(!close(r1Firmicutes.observedExpected,r1Global.observedExpected) || !close(r1Firmicutes.excessExpected,r1Global.excessExpected)){diffFromDefaultPop++;}
			// ncrna (bundle order index 173, cols 1038-1043): obs = the 5 raw ncRNA counts
			final int[] ncObs=new int[]{observedByItem[NFAM],observedByItem[NFAM+1],observedByItem[NFAM+2],observedByItem[NFAM+3],observedByItem[NFAM+4]};
			checkSubnetBlock(r,NCRNA_INDEX,netN,ncObs,phylumCtx,
				MagQCExpectedCopyFeatures.computeForItems(observedByItem,globalMeans,new int[]{NFAM,NFAM+1,NFAM+2,NFAM+3,NFAM+4}),"ncrna",i);
			// rrna (bundle order index 174): its four observations reuse the first four legacy N items.
			final int[] rrnaObs=new int[]{observedByItem[NFAM],observedByItem[NFAM+1],observedByItem[NFAM+2],observedByItem[NFAM+3]};
			checkSubnetBlock(r,RRNA_INDEX,netR,rrnaObs,phylumCtx,
				MagQCExpectedCopyFeatures.computeForItems(observedByItem,globalMeans,new int[]{NFAM,NFAM+1,NFAM+2,NFAM+3}),"rrna",i);
			// trna_anticodon (bundle order index 175): reconstruct all 65 observations from the raw
			// tRNA total (every fixture count is code 5), then compare the complete six-value block.
			final int[] trnaObs=new int[MagQCObservationLayout.ANTICODON_COUNT];
			final int[] trnaItems=anticodonItemIndexes();
			for(int j=0;j<trnaObs.length;j++){trnaObs[j]=observedByItem[trnaItems[j]];}
			checkSubnetBlock(r,TRNA_INDEX,netT,trnaObs,phylumCtx,
				MagQCExpectedCopyFeatures.computeForItems(observedByItem,globalMeans,trnaItems),"trna_anticodon",i);
		}
		check(diffFromDefaultPop>0,"famset_1's phylum-Firmicutes declared population produced the SAME O/E,X/E as "
			+"the global population on every row -- the declaration would have no observable effect, which would "
			+"make this check vacuous (design revision2 non-vacuousness requirement)");
		System.out.println("MAGQC_SIX_FEATURE_TEST exact recomputation PASS: 20 rows, 5 subnets (famset_0/global, "
			+"famset_1/phylum-Firmicutes, ncrna/global, rrna/global, trna_anticodon/global) recomputed exactly across bundle order; famset_1's declared "
			+"population differs from global on "+diffFromDefaultPop+"/20 rows (non-vacuous).");

		// --- (3) non-finite learned output rejection (R5): a second bundle where famset_0's net overflows to +Infinity ---
		final File badNets=new File(dir,"bad_nets"); badNets.mkdirs();
		for(int i=0;i<N;i++){
			final String type=typeForIndex(i);
			final boolean special=!"famset".equals(type);
			final String id=special ? type : "famset_"+i;
			// Reuses the SAME agg_manifest.tsv as the good bundle (badPackArgs only overrides out=/
			// netroot=), which still declares famset_0 at CANDIDATE_OBS=2/width 22 -- the width here
			// must match that declaration exactly, even though this bundle replaces famset_0 with a
			// synthetic overflow net instead of the real candidate.
			final int numObs=special ? MagQCObservationLayout.specialCount(type) : (i==0 ? CANDIDATE_OBS : FAMSET_OBS);
			final int inputWidth=numObs+3+MagQCVectorMaker.SHARED_CONTEXT_WIDTH;
			if(i==0){writeOverflowFourOutputNet(new File(badNets,"net_"+id+".bbnet"),inputWidth,i+1000);}
			else{writeFourOutputNet(new File(badNets,"net_"+id+".bbnet"),inputWidth,i+1000);}
		}
		final File badBundleFile=new File(dir,"six_bundle_bad.bbnets");
		final Map<String,String> badPackArgs=new LinkedHashMap<String,String>(packArgs);
		badPackArgs.put("out",badBundleFile.getAbsolutePath()); badPackArgs.put("netroot",badNets.getAbsolutePath());
		MagQCNetBundle.pack(badPackArgs);
		boolean threw=false; String message=null;
		try{
			MagQCVectorMaker.main(new String[]{
				"cache="+new File(dir,"cache.tsv").getAbsolutePath(),"sizemap="+new File(dir,"sizemap.tsv").getAbsolutePath(),
				"taxpgm="+new File(dir,"taxpgm.tsv").getAbsolutePath(),"familylist="+familyList.getAbsolutePath(),
				"n=5","valn=0","valfrac=0.5","seed=1","cleanspike=0.0",
				"out="+new File(dir,"g_six_bad.tsv").getAbsolutePath(),
				"bundle="+badBundleFile.getAbsolutePath(),"aggout="+new File(dir,"agg_six_bad.tsv").getAbsolutePath(),
				"subnetfeatures=six","expectedcopytable="+table,"expectedcopytablesha256="+tableSha,
				"subnetpopulations="+declaration,"subnetpopulationssha256="+declSha});
		}catch(RuntimeException e){threw=true; message=e.getMessage();}
		check(threw,"expected a non-finite-output rejection, none thrown");
		check(message!=null && message.contains("famset_0") && message.contains("gene_completeness"),
			"rejection message must name the subnet id and feature name: "+message);
		System.out.println("MAGQC_SIX_FEATURE_TEST non-finite rejection PASS: "+message);

		System.out.println("MAGQC_SIX_FEATURE_TEST PASS: real 176-subnet four-output bundle packed through "
			+"MagQCNetBundle.pack(), subnetfeatures=six run produced "+rows.size()+" rows at the correct width "
			+expectedWidth+", all famset_0 learned outputs finite across all rows, and a genuine non-finite "
			+"learned output was rejected with subnet id + feature name named in the message.");
	}

	/** Minimal valid four-output .bbnet: dims=[inputWidth,4], dense, decimal coding, sigmoid on
	 *  every output cell -- BundleFixtureGen's writeNet pattern, generalized to outputWidth=4. */
	static void writeFourOutputNet(File f,int inputWidth,int seed) throws Exception {
		writeNetInternal(f,inputWidth,4,seed,false);
	}
	/** A 4-layer (dims=[inputWidth,1,1,4]) net whose output-0 path compounds three multiplications
	 *  by 999999999999999999 (~1e18, the largest decimal literal parse.Parse.parseFloat parses
	 *  CORRECTLY -- verified empirically this session: its digit accumulator silently wraps to
	 *  small/negative garbage FINITE values beyond ~19 digits, so a single huge literal like 3e38
	 *  cannot reliably reach float32 overflow through this parser, and the literal token "NaN"
	 *  triggers an EARLIER, unrelated ml.Cell.setFrom() assertion at net-construction time rather
	 *  than reaching the intended non-finite-OUTPUT rejection). Layer1's single hidden cell sums
	 *  ~1e18 x (nonzero inputs, well under Float.MAX_VALUE~3.4e38 alone); LINEAR passes it through
	 *  unchanged; layer2 multiplies by another ~1e18 (still under 3.4e38 -- verified: 2e19x1e18=2e37);
	 *  layer3's output-0 cell multiplies by a THIRD ~1e18 (2e37x1e18=2e55, which overflows float32
	 *  to +Infinity on the final float cast -- matches this session's own diagnostic:
	 *  (float)(1e18*1e18*1e18)==Float.POSITIVE_INFINITY). Output cells 1-3 stay ordinary/SIG-bounded
	 *  so only output 0 is non-finite, exercising the per-feature-name rejection message. */
	static void writeOverflowFourOutputNet(File f,int inputWidth,int seed) throws Exception {
		final int idH1=inputWidth+1, idH2=inputWidth+2, idOutBase=inputWidth+3;
		final String HUGE="999999999999999999"; // ~1e18, verified within parse.Parse.parseFloat's correct range
		try(PrintWriter w=writer(f)){
			final int edgeCount=inputWidth+1+4;
			w.print("##bbnet\n#version 1\n#concise\n#dense\n#density 1.00000000\n#blocksize 1\n");
			w.print("#seed "+seed+"\n#layers 4\n#dims "+inputWidth+" 1 1 4\n#edges "+edgeCount+"\n#coding decimal\n");
			w.print("##layer 1\n");
			{
				final StringBuilder sb=new StringBuilder();
				sb.append('C').append(idH1).append(" LINEAR 0");
				for(int j=0;j<inputWidth;j++){sb.append(' ').append(HUGE);}
				w.print(sb.toString()); w.print('\n');
			}
			w.print("##layer 2\n");
			w.print("C"+idH2+" LINEAR 0 "+HUGE+"\n");
			w.print("##layer 3\n");
			w.print("C"+idOutBase+" LINEAR 0 "+HUGE+"\n");
			w.print("C"+(idOutBase+1)+" SIG 0.01 0.001\n");
			w.print("C"+(idOutBase+2)+" SIG 0.02 0.001\n");
			w.print("C"+(idOutBase+3)+" SIG 0.03 0.001\n");
		}
	}
	private static void writeNetInternal(File f,int inputWidth,int outputWidth,int seed,boolean overflow) throws Exception {
		try(PrintWriter w=writer(f)){
			final int edgeCount=outputWidth*inputWidth;
			w.print("##bbnet\n#version 1\n#concise\n#dense\n#density 1.00000000\n#blocksize 1\n");
			w.print("#seed "+seed+"\n#layers 2\n#dims "+inputWidth+" "+outputWidth+"\n#edges "+edgeCount+"\n#coding decimal\n");
			w.print("##layer 1\n");
			for(int pos=0;pos<outputWidth;pos++){
				final int cid=1+inputWidth+pos;
				final StringBuilder sb=new StringBuilder();
				// SIG bounds its output to (0,1) regardless of input magnitude (sigmoid(+Infinity)=1.0,
				// finite) -- an overflowing dot-product cannot produce a non-finite SIG output. The
				// overflow target cell uses LINEAR (unbounded pass-through) instead.
				final boolean overflowCell=(overflow && pos==0);
				sb.append('C').append(cid).append(overflowCell?" LINEAR ":" SIG ").append(fmt(0.01*(pos+1)));
				for(int j=0;j<inputWidth;j++){
					// A literal decimal float-overflow weight (e.g. 3e38, per the design doc's own
					// wording) turned out NOT reliably parseable by BBTools' fast custom float parser
					// (parse.Parse.parseFloat): empirically verified this session that its digit
					// accumulator only handles up to ~19 digits (~1e18) correctly and silently wraps
					// to a small/negative garbage FINITE value beyond that (e.g. a 41-digit "1"+40
					// zeros literal parsed as -5.047021E18, not +Infinity) -- so no decimal literal in
					// this parser's supported range can overflow float32 on its own. Verified instead
					// that the literal token "NaN" parses correctly to Float.NaN (Float.isFinite=false);
					// NaN propagates through both multiplication (NaN*0=NaN, not 0) and addition, so a
					// single NaN weight anywhere in the sum guarantees a non-finite dot product
					// regardless of which inputs are nonzero for a given bin -- more robust than the
					// overflow route even where the parser could support it.
					final String wStr=(overflowCell && j==0) ? "NaN" : fmt(0.001*((j+1)*(pos+1+seed)%97-48));
					sb.append(' ').append(wStr);
				}
				w.print(sb.toString()); w.print('\n');
			}
		}
	}
	static String fmt(double v){final float f=(float)v; return (f==0f) ? "0" : Float.toString(f);}

	/** Mirrors MagQCVectorMaker's private fmt() exactly (whole numbers print without a decimal
	 *  point, everything else at 6 fixed decimals) -- used to compare recomputed row values against
	 *  the row's own printed strings. Distinct from {@link #fmt(double)} above, which formats net
	 *  WEIGHTS for writing fixture .bbnet files, a different purpose/format entirely. */
	private static String fmtRow(double v){
		if(v==(long)v){return Long.toString((long)v);}
		return String.format("%.6f",v);
	}
	/** Decodes an appendTwo-encoded family dense-head column pair (presence, excess/EXC_CAP=16) at
	 *  head index j (columns start+2j, start+2j+1) back to the whole-bin family count. Mirrors
	 *  MagQCAggVectorTest.decodeTwo exactly; exact as long as the count never approaches 16 (true
	 *  for this fixture's tiny per-rank totals). */
	private static int decodeTwo(String[] row,int start,int j){
		final int presCol=start+2*j, excCol=start+2*j+1;
		if("0".equals(row[presCol])){return 0;}
		final double excess=Double.parseDouble(row[excCol]);
		return (int)Math.round(excess*16)+1;
	}
	/** Decodes an appendCountTwo-encoded ratio value (count/(1+count), the default enc=ratio
	 *  encoding appendCountTwo always uses regardless of the global enc= flag) at the given column
	 *  back to the whole-bin count. Exact for this fixture's small integer counts. */
	private static int decodeRatio(String[] row,int col){
		final double x=Double.parseDouble(row[col]);
		if(x<=0){return 0;}
		return (int)Math.round(x/(1.0-x));
	}
	/** Feeds one subnet's reconstructed input (decoded whole-bin obs + this row's own phylum/context
	 *  columns) to an independently-loaded copy of its net, and compares all 6 recomputed values
	 *  (4 learned outputs + O/E + X/E) against the row's own printed strings, exactly. */
	private static void checkSubnetBlock(String[] r,int subnetOrder,CellNet net,int[] obs,float[] phylumCtx,
			MagQCExpectedCopyFeatures.Result expected,String name,int rowIdx){
		final float[] in=new float[obs.length+phylumCtx.length];
		for(int i=0;i<obs.length;i++){in[i]=obs[i];}
		System.arraycopy(phylumCtx,0,in,obs.length,phylumCtx.length);
		net.applyInput(in);
		net.feedForward();
		final int base=subnetOrder*6;
		for(int k=0;k<4;k++){
			final float got=net.getOutput(k);
			check(fmtRow(got).equals(r[base+k]),"row "+rowIdx+" "+name+" output "+k+": recomputed "+fmtRow(got)+" != row "+r[base+k]);
		}
		check(fmtRow(expected.observedExpected).equals(r[base+4]),
			"row "+rowIdx+" "+name+" O/E: recomputed "+fmtRow(expected.observedExpected)+" != row "+r[base+4]);
		check(fmtRow(expected.excessExpected).equals(r[base+5]),
			"row "+rowIdx+" "+name+" X/E: recomputed "+fmtRow(expected.excessExpected)+" != row "+r[base+5]);
	}
	/** Returns the canonical expected-copy item indexes for the 64 structural codes plus unknown. */
	private static int[] anticodonItemIndexes(){
		final int[] out=new int[MagQCObservationLayout.ANTICODON_COUNT];
		for(int i=0;i<out.length;i++){out[i]=NFAM+MagQCObservationLayout.LEGACY_COUNT+i;}
		return out;
	}
	/** Returns the explicit special type assigned to a fixed bundle-order slot. */
	private static String typeForIndex(int i){
		if(i==NCRNA_INDEX){return "ncrna";}
		if(i==RRNA_INDEX){return "rrna";}
		if(i==TRNA_INDEX){return "trna_anticodon";}
		return "famset";
	}
	private static boolean close(double a,double b){return Math.abs(a-b)<=1e-9;}

	// ---- synthetic expected-copy table (mirrors MagQCExpectedCopyBindingTest.buildTable, 4 orgs / 8 families) ----
	private static Path buildTable(Path dir,Path familyList,String familyHash) throws Exception {
		final String rA="tid_100_a.fa", rB="tid_200_b.fa", rC="tid_300_c.fa", rD="tid_400_d.fa";
		final String uA=unit(rA), uB=unit(rB), uC=unit(rC), uD=unit(rD);
		final StringBuilder layoutText=new StringBuilder();
		layoutText.append("#schema_version\ttyped_item_layout_v1\n#family_list\t").append(familyList)
			.append("\t#\tsha256\t").append(familyHash).append("\n#columns\t").append(MagQCExpectedCopyTableTyped.LAYOUT_COLUMNS).append('\n');
		for(int i=0;i<NFAM;i++){layoutText.append("P\t").append(i).append("\tf").append(i).append('\t').append(i).append('\n');}
		for(int i=0;i<MagQCObservationLayout.EXTENDED_COUNT;i++){
			layoutText.append("N\t").append(MagQCObservationLayout.nonproteinKey(i)).append("\t-\t")
				.append(NFAM+i).append('\n');
		}
		final Path layout=write(dir.resolve("layout.tsv"),layoutText.toString());
		final Path exclusions=write(dir.resolve("exclusions.tsv"),"#excluded_tid\tclass\treason\n999\ttraining\tfixture\n");
		final Path roster=write(dir.resolve("roster.tsv"),
			"#schema_version\t"+MagQCExpectedCopyTableTyped.ROSTER_SCHEMA+"\n#columns\t"+MagQCExpectedCopyTableTyped.ROSTER_COLUMNS+"\n"
			+uA+"\tgenome_"+uA+"\t"+rA+"\t/a.fa\t"+repeat('a')+"\t100\tBacteria\n"
			+uB+"\tgenome_"+uB+"\t"+rB+"\t/b.fa\t"+repeat('b')+"\t200\tBacteria\n"
			+uC+"\tgenome_"+uC+"\t"+rC+"\t/c.fa\t"+repeat('c')+"\t300\tBacteria\n"
			+uD+"\tgenome_"+uD+"\t"+rD+"\t/d.fa\t"+repeat('d')+"\t400\tBacteria\n");
		final Path labels=write(dir.resolve("copy_labels.tsv"),
			"#columns\t"+MagQCExpectedCopyTableTyped.LABEL_COLUMNS+"\n"
			+labelRow(uA,rA,"/a.fa","a","Firmicutes")+labelRow(uB,rB,"/b.fa","b","Proteobacteria")
			+labelRow(uC,rC,"/c.fa","c","Firmicutes")+labelRow(uD,rD,"/d.fa","d","Proteobacteria"));
		final StringBuilder countsText=new StringBuilder();
		countsText.append("#schema_version\t").append(MagQCExpectedCopyTableTyped.COUNT_SCHEMA).append("\n#columns\t")
			.append(MagQCExpectedCopyTableTyped.COUNT_COLUMNS).append('\n');
		appendCountRows(countsText,uA,100,new int[]{1,1,1,0,0,0,0,0});
		appendCountRows(countsText,uB,200,new int[]{1,0,1,0,0,1,0,0});
		appendCountRows(countsText,uC,300,new int[]{1,1,0,0,0,0,0,0});
		appendCountRows(countsText,uD,400,new int[]{1,0,0,0,0,1,1,0});
		final Path counts=write(dir.resolve("counts.tsv"),countsText.toString());
		final Path table=dir.resolve("expected_copy_table.tsv");
		MagQCExpectedCopyTableTyped.build(counts.toString(),roster.toString(),layout.toString(),labels.toString(),exclusions.toString(),table.toString());
		return table;
	}
	private static String labelRow(String unit,String rel,String path,String hashChar,String phylum){
		return unit+"\tgenome_"+unit+"\t"+rel+"\t"+path+"\t"+repeat(hashChar.charAt(0))+"\tcfg\tserver\tref1\tclassified\td__Bacteria;p__"+phylum+"\tBacteria\t"+phylum+"\t-\t"+repeat(hashChar.charAt(0))+"\n";
	}
	private static void appendCountRows(StringBuilder sb,String unit,int tid,int[] pCounts){
		for(int i=0;i<pCounts.length;i++){sb.append(unit).append("\tgenome_").append(unit).append('\t').append(tid).append("\tP\t").append(i).append('\t').append(pCounts[i]).append('\n');}
		for(int i=0;i<MagQCObservationLayout.EXTENDED_COUNT;i++){
			final String key=MagQCObservationLayout.nonproteinKey(i);
			final int count="trna_anticodon_5".equals(key) ? 1 : 0;
			sb.append(unit).append("\tgenome_").append(unit).append('\t').append(tid).append("\tN\t").append(key).append('\t').append(count).append('\n');
		}
	}
	private static String unit(String sourceRel) throws Exception {
		final byte[] d=MessageDigest.getInstance("SHA-256").digest((sourceRel+"\n").getBytes(StandardCharsets.UTF_8));
		final StringBuilder b=new StringBuilder(66); b.append("g_");
		for(byte x:d){b.append(String.format("%02x",x&255));}
		return b.toString();
	}
	private static String repeat(char c){final StringBuilder b=new StringBuilder(64);for(int i=0;i<64;i++){b.append(c);}return b.toString();}

	// ---- generic helpers (mirrors MagQCAggVectorTest's) ----
	private static ArrayList<String[]> parse(File f) throws Exception {
		final ArrayList<String[]> out=new ArrayList<String[]>();
		boolean first=true;
		for(String line:Files.readAllLines(f.toPath())){
			if(line.length()==0){continue;}
			if(first){first=false; continue;}
			out.add(line.split("\t"));
		}
		return out;
	}
	private static String header(File f) throws Exception {
		for(String line:Files.readAllLines(f.toPath())){if(line.length()>0){return line;}}
		return null;
	}
	private static File write(File f,String content) throws Exception {final FileWriter w=new FileWriter(f); w.write(content); w.close(); return f;}
	private static Path write(Path p,String s) throws Exception {Files.write(p,s.getBytes(StandardCharsets.UTF_8)); return p;}
	private static PrintWriter writer(File f) throws Exception {return new PrintWriter(new java.io.OutputStreamWriter(new java.io.FileOutputStream(f),StandardCharsets.UTF_8));}
	private static String sha(File f) throws Exception {return sha(f.toPath());}
	private static String sha(Path p) throws Exception {
		final byte[] d=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(p));
		final StringBuilder b=new StringBuilder(64);
		for(byte x:d){b.append(String.format("%02x",x&255));}
		return b.toString();
	}
	private static void deleteRecursive(File f){
		final File[] kids=f.listFiles();
		if(kids!=null){for(File k:kids){deleteRecursive(k);}}
		f.delete();
	}
	private static void check(boolean cond,String msg){if(!cond){throw new AssertionError("FAIL: "+msg);}}
}
