package prot;

import java.io.BufferedOutputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import ml.CellNet;
import ml.Function;

/** Focused Java-8 fixture for {@link MagQCExpectedCopyBinding}. All inputs are synthetic. */
public final class MagQCExpectedCopyBindingTest {
	private MagQCExpectedCopyBindingTest(){}

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("expected-copy-binding.");
		try{
			final int numFam=7;
			final Fixture fx=buildTable(dir);
			assertMeans(fx.table_);

			Arrays.fill(Function.TYPE_RATES,0f);
			final byte[] netAlpha=makeNetBytes(new int[]{9,3,1},11L);
			final byte[] netRandom=makeNetBytes(new int[]{4,3,1},13L);
			final byte[] netNcrna=makeNetBytes(new int[]{7,3,1},17L);
			final byte[] netBeta=makeNetBytes(new int[]{4,3,1},19L);
			final int[] dimsAlpha=readDims(netAlpha), dimsRandom=readDims(netRandom);
			final int[] dimsNcrna=readDims(netNcrna), dimsBeta=readDims(netBeta);
			final int seedAlpha=readSeed(netAlpha), seedRandom=readSeed(netRandom);
			final int seedNcrna=readSeed(netNcrna), seedBeta=readSeed(netBeta);

			final SubnetSpec s0=new SubnetSpec("marker_v4b_Alpha","marker_v4b_Alpha","famset",0,7,dimsAlpha[0],1,
					seedAlpha,dimsAlpha,new int[]{0,1,2,3,4,5,6},"marker_v4b_Alpha.bbnet","marker_v4b_Alpha.txt",
					repeat('1'),netAlpha,sha256Hex(netAlpha));
			final SubnetSpec s1=new SubnetSpec("random_v4b_1","random_v4b_1","famset",1,2,dimsRandom[0],1,
					seedRandom,dimsRandom,new int[]{0,3},"random_v4b_1.bbnet","random_v4b_1.txt",
					repeat('2'),netRandom,sha256Hex(netRandom));
			final SubnetSpec s2=new SubnetSpec("ncrna","ncrna","ncrna",2,5,dimsNcrna[0],1,
					seedNcrna,dimsNcrna,null,"ncrna.bbnet","-","-",netNcrna,sha256Hex(netNcrna));
			final SubnetSpec s3good=new SubnetSpec("m90_v4b_Beta","m90_v4b_Beta","famset",3,2,dimsBeta[0],1,
					seedBeta,dimsBeta,new int[]{5,6},"m90_v4b_Beta.bbnet","m90_v4b_Beta.txt",
					repeat('3'),netBeta,sha256Hex(netBeta));
			final SubnetSpec s3badrank=new SubnetSpec("m90_v4b_Beta","m90_v4b_Beta","famset",3,2,dimsBeta[0],1,
					seedBeta,dimsBeta,new int[]{6,7},"m90_v4b_Beta.bbnet","m90_v4b_Beta.txt",
					repeat('3'),netBeta,sha256Hex(netBeta));
			final SubnetSpec[] specsGood={s0,s1,s2,s3good};
			final SubnetSpec[] specsBadRank={s0,s1,s2,s3badrank};

			final String releaseText="marker_v4b_Alpha\tfamset\tsubsets/marker_v4b_Alpha.txt\n"
					+"random_v4b_1\tfamset\tsubsets/random_v4b_1.txt\n"
					+"ncrna\tncrna\t-\n"
					+"m90_v4b_Beta\tfamset\tsubsets/m90_v4b_Beta.txt\n";
			final byte[] releaseBytes=releaseText.getBytes(StandardCharsets.UTF_8);
			final String releaseSha=sha256Hex(releaseBytes);
			final String releaseB64=Base64.getEncoder().encodeToString(releaseBytes);
			final String familyValues="f0,f1,f2,f3,f4,f5,f6";
			final String familySha=sha(fx.family);
			final String altFamilyText="#rank\trep_id\tocc_total\n0\tf0\t1\n1\tf1\t1\n2\tf2\t1\n3\tf3\t1\n4\tf4\t1\n5\tf5\t1\n6\tf6x\t1\n";
			final String altFamilySha=sha256Hex(altFamilyText.getBytes(StandardCharsets.UTF_8));

			final LinkedHashMap<String,String> headerGood=baseHeader(releaseSha,releaseB64,familySha,familyValues);
			final LinkedHashMap<String,String> headerWrongFamily=baseHeader(releaseSha,releaseB64,altFamilySha,familyValues);
			final LinkedHashMap<String,String> headerBadRank=baseHeader(releaseSha,releaseB64,familySha,familyValues);

			final Path bundlePathGood=dir.resolve("bundle.bbnets");
			final Path bundlePathWrongFamily=dir.resolve("bundle_wrongfamily.bbnets");
			final Path bundlePathBadRank=dir.resolve("bundle_badrank.bbnets");
			final MagQCNetBundle bundleGood=buildAndLoad(bundlePathGood,headerGood,specsGood);
			final MagQCNetBundle bundleWrongFamily=buildAndLoad(bundlePathWrongFamily,headerWrongFamily,specsGood);
			final MagQCNetBundle bundleBadRank=buildAndLoad(bundlePathBadRank,headerBadRank,specsBadRank);

			final String declAText="#schema_version\tsubnet_population_declaration_v1\n"
					+"#release_manifest_sha256\t"+releaseSha+"\n"
					+"#columns\tsubnet_id\tpopulation_type\tpopulation_name\n"
					+"marker_v4b_Alpha\tphylum\tAlpha\n"
					+"random_v4b_1\tglobal\t-\n"
					+"ncrna\tglobal\t-\n"
					+"m90_v4b_Beta\tphylum\tAlpha\n";
			final Path declA=write(dir.resolve("declA.tsv"),declAText);
			final String declASha=sha(declA);
			final String declBText=declAText.replace("m90_v4b_Beta\tphylum\tAlpha\n","m90_v4b_Beta\tphylum\tBeta\n");
			final Path declB=write(dir.resolve("declB.tsv"),declBText);
			final String declBSha=sha(declB);

			final String zeros64=repeat('0');
			final Path pWrongRm=write(dir.resolve("decl_wrong_rm.tsv"),declAText.replace(releaseSha,zeros64));
			final String shaWrongRm=sha(pWrongRm);
			final Path pMissingNcrna=write(dir.resolve("decl_missing_ncrna.tsv"),declAText.replace("ncrna\tglobal\t-\n",""));
			final String shaMissingNcrna=sha(pMissingNcrna);
			final Path pExtra=write(dir.resolve("decl_extra.tsv"),declAText+"bogus_subnet\tglobal\t-\n");
			final String shaExtra=sha(pExtra);
			final Path pDupNcrna=write(dir.resolve("decl_dup_ncrna.tsv"),declAText+"ncrna\tglobal\t-\n");
			final String shaDupNcrna=sha(pDupNcrna);
			final Path pGamma=write(dir.resolve("decl_gamma.tsv"),declAText.replace("random_v4b_1\tglobal\t-\n","random_v4b_1\tphylum\tGamma\n"));
			final String shaGamma=sha(pGamma);
			final Path pGlobalX=write(dir.resolve("decl_global_x.tsv"),declAText.replace("random_v4b_1\tglobal\t-\n","random_v4b_1\tglobal\tX\n"));
			final String shaGlobalX=sha(pGlobalX);
			final Path pFoo=write(dir.resolve("decl_foo.tsv"),declAText.replace("#columns\t","#foo\tbar\n#columns\t"));
			final String shaFoo=sha(pFoo);

			// --- 4.4.2 ---
			final MagQCExpectedCopyBinding bA=MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,bundleGood,
					declA.toString(),declASha,numFam);
			check(bA.size()==4,"size==4");
			check(bA.itemWidth()==12,"itemWidth==12");
			check("marker_v4b_Alpha".equals(bA.subnetId(0)),"subnetId(0)");
			check("phylum".equals(bA.populationType(3))&&"Alpha".equals(bA.populationName(3)),"populationType/Name(3)");
			check(bA.familyListSha256().equals(bundleGood.metadata("family_list_sha256")),"familyListSha256");
			check(bA.releaseManifestSha256().equals(bundleGood.metadata("release_manifest_sha256")),"releaseManifestSha256");
			check(bA.tableSha256().equals(fx.tableSha),"tableSha256");

			// --- 4.4.3: worked example ---
			final int[] obs=new int[12];
			obs[0]=1; obs[3]=2; obs[4]=3;
			MagQCExpectedCopyFeatures.Result r=bA.compute(0,obs);
			check(r.expected==7,"worked E");
			check(r.observedCapped==5,"worked O");
			check(r.excess==1,"worked X");
			check(close(r.observedExpected,5.0/7),"worked O/E");
			check(close(r.excessExpected,1.0/7),"worked X/E");

			// --- 4.4.4: fractional (random_v4b_1, global) ---
			MagQCExpectedCopyFeatures.Result rFrac=bA.compute(1,obs);
			check(close(rFrac.expected,4.0/3),"fractional E");
			check(close(rFrac.observedCapped,4.0/3),"fractional O");
			check(close(rFrac.excess,5.0/3),"fractional X");
			check(close(rFrac.observedExpected,1.0),"fractional O/E");
			check(close(rFrac.excessExpected,5.0/4),"fractional X/E");
			final double[] globalMeans=fx.table_.population("global","-").means();
			final MagQCExpectedCopyFeatures.Result crossFrac=MagQCExpectedCopyFeatures.computeForItems(obs,globalMeans,new int[]{0,3});
			check(close(crossFrac.expected,rFrac.expected)&&close(crossFrac.observedCapped,rFrac.observedCapped)
					&&close(crossFrac.excess,rFrac.excess),"fractional cross-check via computeForItems");

			// --- 4.4.5: ncrna on global ---
			obs[7]=2; obs[11]=1;
			final MagQCExpectedCopyFeatures.Result rNcrna=bA.compute(2,obs);
			final MagQCExpectedCopyFeatures.Result crossNcrna=MagQCExpectedCopyFeatures.computeForItems(obs,globalMeans,new int[]{7,8,9,10,11});
			check(close(crossNcrna.expected,rNcrna.expected)&&close(crossNcrna.observedCapped,rNcrna.observedCapped)
					&&close(crossNcrna.excess,rNcrna.excess),"ncrna cross-check via computeForItems");

			// --- 4.4.6: zero-E under declB ---
			final MagQCExpectedCopyBinding bB=MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,bundleGood,
					declB.toString(),declBSha,numFam);
			obs[5]=3;
			final MagQCExpectedCopyFeatures.Result rZero=bB.compute(3,obs);
			check(rZero.expected==0,"zero-E expected");
			check(rZero.excess==3,"zero-E excess");
			check(rZero.observedExpected==0&&rZero.excessExpected==0,"zero-E ratios");

			// --- 4.4.7: explicit-name independence ---
			obs[6]=1;
			final MagQCExpectedCopyFeatures.Result rA3=bA.compute(3,obs);
			final MagQCExpectedCopyFeatures.Result rB3=bB.compute(3,obs);
			check(close(rA3.expected,1.0),"declA subnet3 E==1 (Alpha means 0,1)");
			check(close(rB3.expected,0.0),"declB subnet3 E==0 (Beta means 0,0)");
			check(!close(rA3.expected,rB3.expected),"subnet name must not decide population");
			final double[] alphaMeans=fx.table_.population("phylum","Alpha").means();
			final MagQCExpectedCopyFeatures.Result crossA3=MagQCExpectedCopyFeatures.computeForItems(obs,alphaMeans,new int[]{5,6});
			check(close(crossA3.expected,rA3.expected)&&close(crossA3.observedCapped,rA3.observedCapped)
					&&close(crossA3.excess,rA3.excess),"declA subnet3 cross-check via computeForItems");

			// --- 4.4.8: untracked isolation ---
			obs[0]=99;
			final MagQCExpectedCopyFeatures.Result r99=bA.compute(3,obs);
			obs[0]=0;
			final MagQCExpectedCopyFeatures.Result r0=bA.compute(3,obs);
			check(close(r99.expected,r0.expected)&&close(r99.observedCapped,r0.observedCapped)&&close(r99.excess,r0.excess),
					"subnet3 must ignore untracked item 0");

			// --- 4.4.9: rejections ---
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),zeros64,
					bundleGood,declA.toString(),declASha,numFam);}},"wrong table hash");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,declA.toString(),zeros64,numFam);}},"wrong declaration hash");
			expectFailureContaining(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleWrongFamily,declA.toString(),declASha,numFam);}},"bundle_wrongfamily","family list hash mismatch");
			expectFailureContaining(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleBadRank,declA.toString(),declASha,numFam);}},"bundle_badrank","family rank out of range");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,declA.toString(),declASha,6);}},"numFam=6");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,declA.toString(),declASha,8);}},"numFam=8");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pWrongRm.toString(),shaWrongRm,numFam);}},"declaration wrong release manifest hash");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pMissingNcrna.toString(),shaMissingNcrna,numFam);}},"declaration missing ncrna row");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pExtra.toString(),shaExtra,numFam);}},"declaration extra row");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pDupNcrna.toString(),shaDupNcrna,numFam);}},"declaration duplicate ncrna row");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pGamma.toString(),shaGamma,numFam);}},"declaration unknown population phylum Gamma");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pGlobalX.toString(),shaGlobalX,numFam);}},"declaration global row with name X");
			expectFailure(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,pFoo.toString(),shaFoo,numFam);}},"declaration unknown metadata line");

			// swapped P0/P1 table: open() must still succeed; bind() must fail on canonical order.
			final String tableText=read(fx.table);
			final String[] lines=tableText.split("\n",-1);
			final ArrayList<Integer> dataIdx=new ArrayList<Integer>();
			for(int i=0;i<lines.length;i++){if(lines[i].length()>0&&lines[i].charAt(0)!='#'){dataIdx.add(i);}}
			check(dataIdx.size()==36,"expected 36 data rows (3 populations x 12 items), got "+dataIdx.size());
			for(int b=0;b<3;b++){
				final int i0=dataIdx.get(b*12), i1=dataIdx.get(b*12+1);
				final String tmp=lines[i0]; lines[i0]=lines[i1]; lines[i1]=tmp;
			}
			final Path swapped=write(dir.resolve("table-swapped.tsv"),join(lines));
			final String swappedSha=sha(swapped);
			MagQCExpectedCopyTableTyped.Table swappedTable=null;
			try{
				swappedTable=MagQCExpectedCopyTableTyped.open(swapped.toString(),swappedSha);
			}catch(IllegalArgumentException e){
				throw new RuntimeException("FAIL: open() rejected the swapped table unexpectedly: "+e.getMessage());
			}
			check(swappedTable!=null,"swapped table failed to parse");
			expectFailureContaining(new Runnable(){public void run(){MagQCExpectedCopyBinding.bind(swapped.toString(),swappedSha,
					bundleGood,declA.toString(),declASha,numFam);}},"swapped canonical order","canonical item order mismatch");

			final int[] obsFinal=obs;
			expectFailure(new Runnable(){public void run(){bA.compute(4,obsFinal);}},"compute bad order");
			expectFailure(new Runnable(){public void run(){bA.compute(0,new int[11]);}},"compute width mismatch");

			// --- 4.4.10: immutability ---
			final MagQCExpectedCopyFeatures.Result before=bA.compute(3,obs);
			final double expectedBefore=before.expected, observedBefore=before.observedCapped, excessBefore=before.excess;
			obs[5]=obs[5]+5;
			check(before.expected==expectedBefore&&before.observedCapped==observedBefore&&before.excess==excessBefore,
					"Result fields must not change after external array mutation");
			final MagQCExpectedCopyFeatures.Result after=bA.compute(3,obs);
			check(after.excess!=excessBefore||after.observedCapped!=observedBefore,
					"fresh compute must reflect obs mutation - stale caching suspected");

			// Declaration file order cannot change bundle order or computed features.
			final String[] declLines=declAText.trim().split("\n");
			final String shuffledText=declLines[0]+"\n"+declLines[1]+"\n"+declLines[2]+"\n"
					+declLines[6]+"\n"+declLines[5]+"\n"+declLines[4]+"\n"+declLines[3]+"\n";
			final Path shuffledDecl=write(dir.resolve("decl-shuffled.tsv"),shuffledText);
			final String shuffledHash=sha(shuffledDecl);
			final MagQCExpectedCopyBinding shuffledBinding=MagQCExpectedCopyBinding.bind(fx.table.toString(),fx.tableSha,
					bundleGood,shuffledDecl.toString(),shuffledHash,numFam);
			check(shuffledHash.equals(shuffledBinding.declarationSha256()),"verified declaration identity retained");
			for(int i=0;i<bA.size();i++){
				check(bA.subnetId(i).equals(shuffledBinding.subnetId(i)),"declaration order must not change bundle order");
				final MagQCExpectedCopyFeatures.Result original=bA.compute(i,obs), shuffled=shuffledBinding.compute(i,obs);
				check(close(original.expected,shuffled.expected)&&close(original.observedCapped,shuffled.observedCapped)
						&&close(original.excess,shuffled.excess),"shuffled declaration feature parity");
			}

			// Bundle exposes its rank array; the binding must own an independent snapshot.
			final MagQCExpectedCopyFeatures.Result snapshotBefore=bA.compute(3,obs);
			final int originalRank=bundleGood.subnet(3).familyRanks[0];
			bundleGood.subnet(3).familyRanks[0]=0;
			final MagQCExpectedCopyFeatures.Result snapshotAfter=bA.compute(3,obs);
			check(close(snapshotBefore.expected,snapshotAfter.expected)
					&&close(snapshotBefore.observedCapped,snapshotAfter.observedCapped)
					&&close(snapshotBefore.excess,snapshotAfter.excess),"binding must snapshot bundle ranks");
			bundleGood.subnet(3).familyRanks[0]=originalRank;

			// --- 4.4.1: open survives source deletion; load fails ---
			Files.deleteIfExists(fx.counts);
			Files.deleteIfExists(fx.roster);
			Files.deleteIfExists(fx.layout);
			Files.deleteIfExists(fx.labels);
			Files.deleteIfExists(fx.exclusions);
			Files.deleteIfExists(fx.family);
			final MagQCExpectedCopyTableTyped.Table reopened=MagQCExpectedCopyTableTyped.open(fx.table.toString(),fx.tableSha);
			check(reopened!=null,"open must survive source deletion");
			boolean loadFailed=false;
			try{
				MagQCExpectedCopyTableTyped.load(fx.table.toString(),fx.counts.toString(),fx.roster.toString(),
						fx.layout.toString(),fx.labels.toString(),fx.exclusions.toString());
			}catch(RuntimeException e){
				loadFailed=true;
			}
			check(loadFailed,"load() must fail once sources are deleted");

			System.out.println("PASS MagQCExpectedCopyBindingTest");
		}finally{deleteTree(dir);}
	}

	private static LinkedHashMap<String,String> baseHeader(String releaseSha,String releaseB64,String familySha,String familyValues){
		final LinkedHashMap<String,String> h=new LinkedHashMap<String,String>();
		h.put("schema_version","1");
		h.put("hash_algorithm","SHA-256");
		h.put("subnet_count","4");
		h.put("release_manifest_sha256",releaseSha);
		h.put("release_manifest_values_b64",releaseB64);
		h.put("family_list_sha256",familySha);
		h.put("family_list_values",familyValues);
		return h;
	}

	private static MagQCNetBundle buildAndLoad(Path out,LinkedHashMap<String,String> header,SubnetSpec[] specs) throws Exception{
		writeBundleFile(out,header,specs);
		return MagQCNetBundle.load(out);
	}

	private static void writeBundleFile(Path out,LinkedHashMap<String,String> header,SubnetSpec[] specs) throws Exception{
		final String canonical=canonicalHash(header,specs);
		final LinkedHashMap<String,String> finalHeader=new LinkedHashMap<String,String>(header);
		finalHeader.put("canonical_payload_sha256",canonical);
		final StringBuilder sb=new StringBuilder();
		sb.append("MAGQC_BBNETS_V1\n");
		for(Map.Entry<String,String> e:finalHeader.entrySet()){
			if("release_manifest_values_b64".equals(e.getKey())){writeWrapped(sb,e.getKey(),e.getValue());}
			else{sb.append(e.getKey()).append('=').append(e.getValue()).append('\n');}
		}
		for(SubnetSpec s:specs){
			sb.append("##subnet\n");
			field(sb,"id",s.id); field(sb,"name",s.name); field(sb,"type",s.type); field(sb,"order",s.order);
			field(sb,"num_obs",s.numObs); field(sb,"expected_inputs",s.expectedInputs); field(sb,"expected_outputs",s.expectedOutputs);
			field(sb,"seed",s.seed); field(sb,"dims",joinInts(s.dims,","));
			field(sb,"family_ranks",ranksToken(s.familyRanks));
			field(sb,"ncrna_obs_definition",s.familyRanks==null?"5 ordered ncRNA observations":"-");
			field(sb,"output_units","completeness_or_contamination_regression");
			field(sb,"output_clamp_ratio","[0,2] ratio clamp; one scalar output");
			field(sb,"loose_net_name",s.looseNetName); field(sb,"loose_subset_name",s.looseSubsetName);
			field(sb,"subset_sha256",s.subsetSha256); field(sb,"net_sha256",s.netSha256); field(sb,"net_bytes",s.netBytes.length);
			writeWrapped(sb,"net_b64",Base64.getEncoder().encodeToString(s.netBytes));
			sb.append("##endsubnet\n");
		}
		final OutputStream fos=new BufferedOutputStream(Files.newOutputStream(out));
		final GZIPOutputStream gz=new GZIPOutputStream(fos);
		gz.write(sb.toString().getBytes(StandardCharsets.UTF_8));
		gz.close();
	}

	/** Mirrors MagQCNetBundle.canonicalHash; if it drifts, load() rejects this fixture and the test fails loud, which is intended. */
	private static String canonicalHash(LinkedHashMap<String,String> header,SubnetSpec[] specs){
		final StringBuilder b=new StringBuilder();
		b.append("MAGQC_BBNETS_CANONICAL_V1\n");
		for(Map.Entry<String,String> e:header.entrySet()){
			if("canonical_payload_sha256".equals(e.getKey())||"input_root_recorded".equals(e.getKey())
					||"subset_root_recorded".equals(e.getKey())){continue;}
			b.append(e.getKey()).append('=').append(e.getValue()).append('\n');
		}
		for(SubnetSpec s:specs){
			b.append("subnet\t").append(s.order).append('\t').append(s.id).append('\t').append(s.name).append('\t').append(s.type).append('\t')
				.append(s.numObs).append('\t').append(s.expectedInputs).append('\t').append(s.expectedOutputs).append('\t').append(s.seed).append('\t')
				.append(joinInts(s.dims,",")).append('\t').append(ranksToken(s.familyRanks)).append('\t')
				.append(s.looseNetName).append('\t').append(s.looseSubsetName).append('\t').append(s.subsetSha256).append('\t')
				.append(s.netSha256).append('\t').append(s.netBytes.length).append('\n');
		}
		return sha256Hex(b.toString().getBytes(StandardCharsets.UTF_8));
	}

	private static void field(StringBuilder sb,String k,Object v){sb.append(k).append('=').append(String.valueOf(v)).append('\n');}
	private static void writeWrapped(StringBuilder sb,String key,String value){
		int p=0; boolean first=true;
		while(p<value.length()){
			final int e=Math.min(p+76,value.length());
			sb.append(first?key+"=":"+");
			sb.append(value,p,e);
			sb.append('\n');
			first=false; p=e;
		}
		if(first){sb.append(key).append("=\n");}
	}
	private static String joinInts(int[] x,String sep){final StringBuilder b=new StringBuilder();for(int i=0;i<x.length;i++){if(i>0){b.append(sep);}b.append(x[i]);}return b.toString();}
	private static String ranksToken(int[] x){return x==null?"-":(x.length==0?"empty":joinInts(x,","));}

	private static byte[] makeNetBytes(int[] dims,long seed){
		final CellNet net=new CellNet(dims,seed,1f,0f,1,new ArrayList<String>());
		net.randomize();
		return net.toBytes().toString().getBytes(StandardCharsets.UTF_8);
	}
	private static int[] readDims(byte[] netBytes){
		for(String line:new String(netBytes,StandardCharsets.UTF_8).split("\r?\n")){
			if(line.startsWith("#dims")){
				final String[] parts=line.substring(5).trim().split("\\s+");
				final int[] r=new int[parts.length];
				for(int i=0;i<r.length;i++){r[i]=Integer.parseInt(parts[i]);}
				return r;
			}
		}
		throw new RuntimeException("missing #dims");
	}
	private static int readSeed(byte[] netBytes){
		for(String line:new String(netBytes,StandardCharsets.UTF_8).split("\r?\n")){
			if(line.startsWith("#seed")){return Integer.parseInt(line.substring(5).trim());}
		}
		throw new RuntimeException("missing #seed");
	}

	private static final class SubnetSpec{
		final String id,name,type,looseNetName,looseSubsetName,subsetSha256,netSha256;
		final int order,numObs,expectedInputs,expectedOutputs;
		final long seed;
		final int[] dims,familyRanks;
		final byte[] netBytes;
		SubnetSpec(String id,String name,String type,int order,int numObs,int expectedInputs,int expectedOutputs,
				long seed,int[] dims,int[] familyRanks,String looseNetName,String looseSubsetName,String subsetSha256,
				byte[] netBytes,String netSha256){
			this.id=id; this.name=name; this.type=type; this.order=order; this.numObs=numObs;
			this.expectedInputs=expectedInputs; this.expectedOutputs=expectedOutputs; this.seed=seed;
			this.dims=dims; this.familyRanks=familyRanks; this.looseNetName=looseNetName;
			this.looseSubsetName=looseSubsetName; this.subsetSha256=subsetSha256; this.netBytes=netBytes; this.netSha256=netSha256;
		}
	}

	private static void assertMeans(MagQCExpectedCopyTableTyped.Table t){
		final double[] alphaP={0,1,0,2,3,0,1};
		final double[] alphaN={1,1,0,0,1};
		final double[] betaP={1,0,0,1,0,0,0};
		final double[] betaN={2,0,0,0,0};
		final double[] globalP={1.0/3,1.0/3,0,1,1,0,1.0/3};
		final double[] globalN={1,1.0/3,0,0,1.0/3};
		for(int i=0;i<7;i++){
			check(close(t.expectation("phylum","Alpha","P",Integer.toString(i)),alphaP[i]),"phylum Alpha P"+i);
			check(close(t.expectation("phylum","Beta","P",Integer.toString(i)),betaP[i]),"phylum Beta P"+i);
			check(close(t.expectation("global","-","P",Integer.toString(i)),globalP[i]),"global P"+i);
		}
		for(int i=0;i<5;i++){
			final String key=MagQCExpectedCopyTableTyped.NCRNA[i];
			check(close(t.expectation("phylum","Alpha","N",key),alphaN[i]),"phylum Alpha N "+key);
			check(close(t.expectation("phylum","Beta","N",key),betaN[i]),"phylum Beta N "+key);
			check(close(t.expectation("global","-","N",key),globalN[i]),"global N "+key);
		}
	}

	private static final class Fixture{
		final MagQCExpectedCopyTableTyped.Table table_;
		final String tableSha;
		final Path table,family,counts,roster,layout,labels,exclusions;
		Fixture(MagQCExpectedCopyTableTyped.Table t,String sha,Path table,Path family,Path counts,Path roster,
				Path layout,Path labels,Path exclusions){
			this.table_=t; this.tableSha=sha; this.table=table; this.family=family; this.counts=counts;
			this.roster=roster; this.layout=layout; this.labels=labels; this.exclusions=exclusions;
		}
	}

	private static Fixture buildTable(Path dir) throws Exception{
		final Path family=write(dir.resolve("families.tsv"),
				"#rank\trep_id\tocc_total\n0\tf0\t1\n1\tf1\t1\n2\tf2\t1\n3\tf3\t1\n4\tf4\t1\n5\tf5\t1\n6\tf6\t1\n");
		final String familyHash=sha(family);
		final String rA="tid_101_a.fa", rC="tid_102_c.fa", rD="tid_103_d.fa";
		final String uA=unit(rA), uC=unit(rC), uD=unit(rD);
		final Path layout=write(dir.resolve("layout.tsv"),
				"#schema_version\ttyped_item_layout_v1\n#family_list\t"+family+"\t#\tsha256\t"+familyHash+"\n#columns\t"
						+MagQCExpectedCopyTableTyped.LAYOUT_COLUMNS+"\n"
				+"P\t0\tf0\t0\nP\t1\tf1\t1\nP\t2\tf2\t2\nP\t3\tf3\t3\nP\t4\tf4\t4\nP\t5\tf5\t5\nP\t6\tf6\t6\n"
				+"N\tr16\t-\t7\nN\tr23\t-\t8\nN\tr5\t-\t9\nN\trother\t-\t10\nN\ttrna\t-\t11\n");
		final Path exclusions=write(dir.resolve("exclusions.tsv"),"#excluded_tid\tclass\treason\n999\ttraining\tfixture\n");
		final String rosterHeader="#schema_version\t"+MagQCExpectedCopyTableTyped.ROSTER_SCHEMA+"\n#columns\t"
				+MagQCExpectedCopyTableTyped.ROSTER_COLUMNS+"\n";
		final Path roster=write(dir.resolve("roster.tsv"),rosterHeader
				+uA+"\tgenome_"+uA+"\t"+rA+"\t/assemblies/a.fa\t"+repeat('a')+"\t101\tBacteria\n"
				+uC+"\tgenome_"+uC+"\t"+rC+"\t/assemblies/c.fa\t"+repeat('c')+"\t102\tBacteria\n"
				+uD+"\tgenome_"+uD+"\t"+rD+"\t/assemblies/d.fa\t"+repeat('d')+"\t103\tArchaea\n");
		final String labelHeader="#columns\t"+MagQCExpectedCopyTableTyped.LABEL_COLUMNS+"\n";
		final Path labels=write(dir.resolve("labels.tsv"),labelHeader
				+uA+"\tgenome_"+uA+"\t"+rA+"\t/assemblies/a.fa\t"+repeat('a')+"\tcfg\tserver\tref1\tclassified\td__Bacteria;p__Alpha\tBacteria\tAlpha\t-\t"+repeat('1')+"\n"
				+uC+"\tgenome_"+uC+"\t"+rC+"\t/assemblies/c.fa\t"+repeat('c')+"\tcfg\tserver\tref1\tclassified\td__Bacteria;p__Beta\tBacteria\tBeta\t-\t"+repeat('2')+"\n"
				+uD+"\tgenome_"+uD+"\t"+rD+"\t/assemblies/d.fa\t"+repeat('d')+"\tcfg\tserver\tref1\tpartial\td__Archaea;p__unknown\tArchaea\tunknown\tno_phylum\t"+repeat('3')+"\n"
				+"excluded_999\tgenome_x\ttid_999_excluded.fa\t/assemblies/x.fa\t"+repeat('e')+"\tcfg\tserver\tref1\tunknown\td__unknown;p__unknown\tunknown\tunknown\ttraining_exclusion\t"+repeat('4')+"\n");
		final Path counts=write(dir.resolve("counts.tsv"),"#schema_version\t"+MagQCExpectedCopyTableTyped.COUNT_SCHEMA+"\n#columns\t"
				+MagQCExpectedCopyTableTyped.COUNT_COLUMNS+"\n"
				+uA+"\tgenome_"+uA+"\t101\tP\t1\t1\n"
				+uA+"\tgenome_"+uA+"\t101\tP\t3\t2\n"
				+uA+"\tgenome_"+uA+"\t101\tP\t4\t3\n"
				+uA+"\tgenome_"+uA+"\t101\tP\t6\t1\n"
				+uA+"\tgenome_"+uA+"\t101\tN\tr16\t1\n"
				+uA+"\tgenome_"+uA+"\t101\tN\tr23\t1\n"
				+uA+"\tgenome_"+uA+"\t101\tN\tr5\t0\n"
				+uA+"\tgenome_"+uA+"\t101\tN\trother\t0\n"
				+uA+"\tgenome_"+uA+"\t101\tN\ttrna\t1\n"
				+uC+"\tgenome_"+uC+"\t102\tP\t0\t1\n"
				+uC+"\tgenome_"+uC+"\t102\tP\t3\t1\n"
				+uC+"\tgenome_"+uC+"\t102\tN\tr16\t2\n"
				+uC+"\tgenome_"+uC+"\t102\tN\tr23\t0\n"
				+uC+"\tgenome_"+uC+"\t102\tN\tr5\t0\n"
				+uC+"\tgenome_"+uC+"\t102\tN\trother\t0\n"
				+uC+"\tgenome_"+uC+"\t102\tN\ttrna\t0\n"
				+uD+"\tgenome_"+uD+"\t103\tN\tr16\t0\n"
				+uD+"\tgenome_"+uD+"\t103\tN\tr23\t0\n"
				+uD+"\tgenome_"+uD+"\t103\tN\tr5\t0\n"
				+uD+"\tgenome_"+uD+"\t103\tN\trother\t0\n"
				+uD+"\tgenome_"+uD+"\t103\tN\ttrna\t0\n");
		final Path table=dir.resolve("table.tsv");
		final MagQCExpectedCopyTableTyped.Table t=MagQCExpectedCopyTableTyped.build(counts.toString(),roster.toString(),
				layout.toString(),labels.toString(),exclusions.toString(),table.toString());
		final String tableSha=sha(table);
		return new Fixture(t,tableSha,table,family,counts,roster,layout,labels,exclusions);
	}

	private static String read(Path p) throws Exception{return new String(Files.readAllBytes(p),StandardCharsets.UTF_8);}
	private static Path write(Path p,String s) throws Exception{Files.write(p,s.getBytes(StandardCharsets.UTF_8));return p;}
	private static String unit(String sourceRel) throws Exception{
		final byte[] d=MessageDigest.getInstance("SHA-256").digest((sourceRel+"\n").getBytes(StandardCharsets.UTF_8));
		final StringBuilder b=new StringBuilder(66); b.append("g_");
		for(byte x:d){b.append(String.format("%02x",x&255));}
		return b.toString();
	}
	private static String repeat(char c){final StringBuilder b=new StringBuilder(64);for(int i=0;i<64;i++){b.append(c);}return b.toString();}
	private static String sha(Path p) throws Exception{return sha256Hex(Files.readAllBytes(p));}
	private static String sha256Hex(byte[] bytes){
		try{
			final byte[] d=MessageDigest.getInstance("SHA-256").digest(bytes);
			final StringBuilder b=new StringBuilder(64);
			for(byte x:d){b.append(String.format("%02x",x&255));}
			return b.toString();
		}catch(Exception e){throw new RuntimeException(e);}
	}
	private static String join(String[] lines){final StringBuilder b=new StringBuilder();for(int i=0;i<lines.length;i++){if(i>0){b.append('\n');}b.append(lines[i]);}return b.toString();}
	private static boolean close(double a,double b){return Math.abs(a-b)<=1e-12;}
	private static void check(boolean ok,String message){if(!ok){throw new RuntimeException("FAIL: "+message);}}
	private static void expectFailure(Runnable r,String what){try{r.run();throw new RuntimeException("FAIL: did not reject "+what);}catch(IllegalArgumentException expected){}}
	private static void expectFailureContaining(Runnable r,String what,String fragment){
		try{r.run();throw new RuntimeException("FAIL: did not reject "+what);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null&&expected.getMessage().indexOf(fragment)>=0,
					"specific failure for "+what+": "+expected.getMessage());
		}
	}
	private static void deleteTree(Path p) throws Exception{
		if(!Files.exists(p)){return;}
		final java.nio.file.DirectoryStream<Path> ds=Files.newDirectoryStream(p);
		try{for(Path child:ds){if(Files.isDirectory(child)){deleteTree(child);}else{Files.deleteIfExists(child);}}}
		finally{ds.close();}
		Files.deleteIfExists(p);
	}
}
