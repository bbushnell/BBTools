package prot;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.util.*;
import java.util.zip.*;

/** Focused contract test for the schema_version=2 (opt-in four-output subnet) increment.
 * Reuses a mixed (175 legacy + 1 four-output) bundle already packed by the real pack() path;
 * this test never weakens a rejection to make a hand-authored fixture pass -- every expectBad
 * case here is a genuine, real defect the loader must catch. */
// Packaged by Yoimiya from the accepted four-output increment fixture.
// Usage: java -ea prot.MagQCNetBundleFourOutputTest MIXED_BUNDLE
public final class MagQCNetBundleFourOutputTest {
	static String FOUR_ID="famset_0", LEGACY_ID="famset_1";

	public static void main(String[] args) throws Exception {
		if(args.length!=1){throw new IllegalArgumentException("Usage: java -ea prot.MagQCNetBundleFourOutputTest MIXED_BUNDLE");}
		String mixedPath=args[0];
		byte[] plain=gunzip(Files.readAllBytes(Paths.get(mixedPath)));

		testDefaultLoadRejectsSchema2(plain);
		testLoadMultiOutputAcceptsAndScores(mixedPath);
		testScoreThrowsForFourOutput(mixedPath);
		testExactNamesUnitsContract(plain);
		testOutputSemanticsTampering(plain);
		testOutputClampRatioTampering(plain);
		testSchemaEntryCountDisagreement(plain);
		testTamperedButHashConsistent(plain);
		testConcurrentScoreAll(mixedPath);
		testRelocationInvariance(mixedPath);
		System.out.println("MAGQC_FOUR_OUTPUT_CONTRACT PASS");
	}

	static void testDefaultLoadRejectsSchema2(byte[] plain) throws Exception {
		try { MagQCNetBundle.load(writeTemp(gzip(plain))); throw new AssertionError("load() accepted schema2"); }
		catch(IOException expected) { check(expected.getMessage().contains("loadMultiOutput"),"rejection names loadMultiOutput"); }
	}

	static void testLoadMultiOutputAcceptsAndScores(String mixedPath) throws Exception {
		MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(mixedPath));
		check(b.size()==176,"mixed bundle size");
		check("2".equals(b.metadata("schema_version")),"schema_version=2");
		MagQCNetBundle.Subnet four=b.subnet(FOUR_ID);
		check(four.expectedOutputs==4,"four-output entry declares 4");
		check(Arrays.equals(four.outputNames(),new String[]{"gene_completeness","gene_contamination",
				"abs_residual_gene_completeness","abs_residual_gene_contamination"}),"exact frozen names");
		float[] in=vec(four.expectedInputs,1);
		float[] out=four.scoreAll(in);
		check(out.length==4,"scoreAll returns 4 values");
		float[] out2=four.scoreAll(in);
		check(Arrays.equals(out,out2),"scoreAll repeatable");
		out2[0]=-999f; // mutate the second snapshot
		check(out[0]!=-999f,"scoreAll returns an INDEPENDENT snapshot, not a shared/aliased array");
	}

	static void testScoreThrowsForFourOutput(String mixedPath) throws Exception {
		MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(mixedPath));
		MagQCNetBundle.Subnet four=b.subnet(FOUR_ID);
		try { four.score(vec(four.expectedInputs,1)); throw new AssertionError("score() did not throw"); }
		catch(IllegalStateException expected) { }
		// A legacy entry in the SAME mixed bundle is unaffected.
		MagQCNetBundle.Subnet legacy=b.subnet(LEGACY_ID);
		check(legacy.expectedOutputs==1,"legacy entry still declares 1");
		float s=legacy.score(vec(legacy.expectedInputs,2));
		float[] all=legacy.scoreAll(vec(legacy.expectedInputs,2));
		check(all.length==1 && all[0]==s,"legacy score()/scoreAll() agree");
	}

	static void testExactNamesUnitsContract(byte[] plain) throws Exception {
		String text=new String(plain,StandardCharsets.UTF_8);
		String origNamesB64=fieldValue(text,"output_names");
		String decoded=new String(Base64.getDecoder().decode(origNamesB64),StandardCharsets.UTF_8);
		// wrong single name (typo in the first entry)
		expectRejectAfterRepack(text,decoded,decoded.replaceFirst("gene_completeness","gene_completenes"),"typo'd name");
		// wrong order (swap entries 0 and 1)
		String[] parts=decoded.split("\n",-1);
		String swapped=parts[1]+"\n"+parts[0]+"\n"+parts[2]+"\n"+parts[3];
		expectRejectAfterRepack(text,decoded,swapped,"reordered names");
		// missing entry (3 instead of 4)
		String missing=parts[0]+"\n"+parts[1]+"\n"+parts[2];
		expectRejectAfterRepack(text,decoded,missing,"missing name entry");
		// extra entry (5 instead of 4)
		String extra=decoded+"\nextra_bogus_name";
		expectRejectAfterRepack(text,decoded,extra,"extra name entry");
		// duplicate entry (4 entries, but two identical, so it can't equal the frozen 4-tuple)
		String duplicate=parts[0]+"\n"+parts[0]+"\n"+parts[2]+"\n"+parts[3];
		expectRejectAfterRepack(text,decoded,duplicate,"duplicate name entry");
	}

	/** Replaces the output_names field's content, recomputes the canonical_payload_sha256 so the
	 * file is INTERNALLY self-consistent (a naive hash-only check would pass), and confirms
	 * validate() still rejects it because the decoded names don't match the frozen external
	 * contract. This is the adversarial case root specifically asked for: hash consistency alone
	 * must not be mistaken for content correctness. */
	static void expectRejectAfterRepack(String origText,String origDecoded,String newDecoded,String label) throws Exception {
		String origB64=Base64.getEncoder().encodeToString(origDecoded.getBytes(StandardCharsets.UTF_8));
		String newB64=Base64.getEncoder().encodeToString(newDecoded.getBytes(StandardCharsets.UTF_8));
		String tampered=origText.replaceFirst("output_names="+java.util.regex.Pattern.quote(origB64),"output_names="+newB64);
		check(!tampered.equals(origText),label+": replacement actually matched something");
		String rehashed=recomputeCanonicalHash(tampered);
		try {
			MagQCNetBundle.loadMultiOutput(writeTemp(gzip(rehashed.getBytes(StandardCharsets.UTF_8))));
			throw new AssertionError(label+": was accepted despite tampered output_names (even with a self-consistent recomputed hash)");
		} catch(IOException expected) {
			check(expected.getMessage()!=null && expected.getMessage().contains("output_names"),
					label+": rejected for the wrong reason: "+expected.getMessage());
		}
	}

	/** Bundle-level output_semantics: a forger changes the string AND recomputes
	 * canonical_payload_sha256 (internally self-consistent) -- validate() must still reject
	 * because the value doesn't equal the frozen SCHEMA2_OUTPUT_SEMANTICS literal, independent
	 * of hash self-consistency (root's review, 2026-09-09). */
	static void testOutputSemanticsTampering(byte[] plain) throws Exception {
		String text=new String(plain,StandardCharsets.UTF_8);
		String origLine=grepLine(text,"output_semantics=");
		String tampered=text.replaceFirst(java.util.regex.Pattern.quote(origLine),
				"output_semantics=plausible but wrong text");
		check(!tampered.equals(text),"output_semantics replacement actually matched");
		String rehashed=recomputeCanonicalHash(tampered);
		try { MagQCNetBundle.loadMultiOutput(writeTemp(gzip(rehashed.getBytes(StandardCharsets.UTF_8))));
			throw new AssertionError("tampered output_semantics (with self-consistent recomputed hash) was accepted"); }
		catch(IOException expected) { check(expected.getMessage().contains("output_semantics"),"rejected for the wrong reason: "+expected.getMessage()); }
	}

	/** Per-subnet output_clamp_ratio on a four-output entry: same self-consistent-forgery shape. */
	static void testOutputClampRatioTampering(byte[] plain) throws Exception {
		String text=new String(plain,StandardCharsets.UTF_8);
		String origLine=grepLine(text,"output_clamp_ratio=");
		String tampered=text.replaceFirst(java.util.regex.Pattern.quote(origLine),
				"output_clamp_ratio=[0,2] ratio clamp; one scalar output");
		check(!tampered.equals(text),"output_clamp_ratio replacement actually matched");
		String rehashed=recomputeCanonicalHash(tampered);
		try { MagQCNetBundle.loadMultiOutput(writeTemp(gzip(rehashed.getBytes(StandardCharsets.UTF_8))));
			throw new AssertionError("tampered output_clamp_ratio (with self-consistent recomputed hash) was accepted"); }
		catch(IOException expected) { check(expected.getMessage().contains("output_clamp_ratio"),"rejected for the wrong reason: "+expected.getMessage()); }
	}
	static String grepLine(String text,String prefix) {
		for(String l:text.split("\n",-1)) if(l.startsWith(prefix)) return l;
		throw new RuntimeException("field not found: "+prefix);
	}

	static void testSchemaEntryCountDisagreement(byte[] plain) throws Exception {
		String text=new String(plain,StandardCharsets.UTF_8);
		// schema_version=1 but a block still declares expected_outputs=4 (and carries names) --
		// simulate by forcing schema_version back to 1 without touching the four-output block.
		String forced=text.replaceFirst("schema_version=2","schema_version=1");
		String rehashed=recomputeCanonicalHash(forced);
		try { MagQCNetBundle.loadMultiOutput(writeTemp(gzip(rehashed.getBytes(StandardCharsets.UTF_8))));
			throw new AssertionError("schema_version=1 with a four-output block was accepted"); }
		catch(IOException expected) { check(expected.getMessage().contains("four-output"),"schema/entry disagreement message"); }
	}

	static void testTamperedButHashConsistent(byte[] plain) throws Exception {
		// Already covered by expectRejectAfterRepack above (every case there recomputes the
		// canonical hash before asserting rejection) -- this method exists as an explicit,
		// separately-named checkpoint so the coverage list is traceable one-to-one with root's
		// requirements, not folded silently into another test's side effect.
		System.out.println("testTamperedButHashConsistent: covered by testExactNamesUnitsContract's expectRejectAfterRepack (every case recomputes canonical_payload_sha256 before asserting rejection)");
	}

	static void testConcurrentScoreAll(String mixedPath) throws Exception {
		MagQCNetBundle b=MagQCNetBundle.loadMultiOutput(Paths.get(mixedPath));
		MagQCNetBundle.Subnet four=b.subnet(FOUR_ID);
		float[] in=vec(four.expectedInputs,3);
		float[] serial=four.scoreAll(in);
		java.util.concurrent.ExecutorService pool=java.util.concurrent.Executors.newFixedThreadPool(4);
		ArrayList<java.util.concurrent.Future<float[]>> fs=new ArrayList<>();
		for(int i=0;i<16;i++) fs.add(pool.submit(() -> four.scoreAll(in)));
		pool.shutdown();
		for(java.util.concurrent.Future<float[]> f:fs) check(Arrays.equals(serial,f.get()),"concurrent scoreAll matches serial");
	}

	static void testRelocationInvariance(String mixedPath) throws Exception {
		MagQCNetBundle orig=MagQCNetBundle.loadMultiOutput(Paths.get(mixedPath));
		// LIMITATION (root's review, 2026-09-09): this only copies an ALREADY-PACKED file's bytes
		// to a new path and re-hashes -- it does NOT exercise the actual
		// input_root_recorded/subset_root_recorded exclusion mechanism, which requires re-PACKING
		// from a genuinely different netroot/subsetroot and confirming the canonical hash is
		// unchanged. This is a weaker check than that (trivially true: an unmodified file's own
		// stored hash does not change when the file is copied). Root is running the real
		// relocated-repack check independently; this test is kept only as a cheap sanity check
		// that loadMultiOutput() itself is path-independent, not as relocation-invariance evidence.
		Path moved=Files.createTempFile("magqc-mixed-relocated-",".bbnets");
		Files.write(moved,Files.readAllBytes(Paths.get(mixedPath)));
		MagQCNetBundle reloaded=MagQCNetBundle.loadMultiOutput(moved);
		check(orig.metadata("canonical_payload_sha256").equals(reloaded.metadata("canonical_payload_sha256")),
				"relocation does not change canonical hash");
		Files.deleteIfExists(moved);
	}

	// ---- helpers ----
	static float[] vec(int n,int seed){float[] v=new float[n];for(int i=0;i<n;i++)v[i]=(float)((i*13+seed*7)%97)/96f;return v;}
	static String fieldValue(String text,String key){
		int i=text.indexOf(key+"=");
		if(i<0) throw new RuntimeException("field not found: "+key);
		int end=text.indexOf('\n',i);
		return text.substring(i+key.length()+1,end);
	}
	static Path writeTemp(byte[] gz) throws IOException {
		Path p=Files.createTempFile("magqc-fourout-test-",".bbnets");
		Files.write(p,gz);
		return p;
	}
	static byte[] gunzip(byte[] x) throws Exception {
		GZIPInputStream in=new GZIPInputStream(new ByteArrayInputStream(x)); ByteArrayOutputStream b=new ByteArrayOutputStream(); byte[] q=new byte[8192]; int n;
		while((n=in.read(q))>=0) if(n>0)b.write(q,0,n); in.close(); return b.toByteArray();
	}
	static byte[] gzip(byte[] x) throws Exception { ByteArrayOutputStream b=new ByteArrayOutputStream(); GZIPOutputStream out=new GZIPOutputStream(b); out.write(x); out.close(); return b.toByteArray(); }

	/** Recomputes canonical_payload_sha256 over the tampered text using the SAME algorithm
	 * MagQCNetBundle.canonicalHash() uses (V2 domain, metadata lines minus the two excluded
	 * provenance fields and canonical_payload_sha256 itself, then one tab-joined line per
	 * subnet block including output_names/output_units for four-output entries) -- an
	 * independent reimplementation, not a call into the class under test, so this genuinely
	 * proves content-contract validation is separate from hash self-consistency. */
	static String recomputeCanonicalHash(String text) throws Exception {
		String[] lines=text.split("\n",-1);
		LinkedHashMap<String,String> md=new LinkedHashMap<>();
		int p=1; // skip magic
		while(p<lines.length && !"##subnet".equals(lines[p])) {
			String l=lines[p++];
			if(l.length()==0||l.startsWith("#")) continue;
			int eq=l.indexOf('=');
			if(eq<=0||l.startsWith("+")) continue;
			String k=l.substring(0,eq), v=l.substring(eq+1);
			if("release_manifest_values_b64".equals(k)) { StringBuilder b=new StringBuilder(v); while(p<lines.length && lines[p].startsWith("+")) b.append(lines[p++].substring(1)); v=b.toString(); }
			md.put(k,v);
		}
		StringBuilder canon=new StringBuilder();
		canon.append("2".equals(md.get("schema_version")) ? "MAGQC_BBNETS_CANONICAL_V2\n" : "MAGQC_BBNETS_CANONICAL_V1\n");
		for(Map.Entry<String,String> e:md.entrySet()) {
			if("canonical_payload_sha256".equals(e.getKey())||"input_root_recorded".equals(e.getKey())||"subset_root_recorded".equals(e.getKey())) continue;
			canon.append(e.getKey()).append('=').append(e.getValue()).append('\n');
		}
		while(p<lines.length) {
			while(p<lines.length && (lines[p].length()==0 || (lines[p].startsWith("#") && !"##subnet".equals(lines[p])))) p++;
			if(p>=lines.length) break;
			if(!"##subnet".equals(lines[p++])) throw new RuntimeException("expected ##subnet");
			LinkedHashMap<String,String> f=new LinkedHashMap<>();
			while(p<lines.length && !"##endsubnet".equals(lines[p])) {
				String l=lines[p++]; if(l.length()==0||l.startsWith("#")) continue;
				int eq=l.indexOf('='); if(eq<=0) continue;
				String k=l.substring(0,eq), v=l.substring(eq+1);
				if("net_b64".equals(k)||"output_names".equals(k)||"output_units".equals(k)||"release_manifest_values_b64".equals(k)) {
					StringBuilder b=new StringBuilder(v); while(p<lines.length && lines[p].startsWith("+")) b.append(lines[p++].substring(1)); v=b.toString();
				}
				f.put(k,v);
			}
			p++; // ##endsubnet
			int expOut=Integer.parseInt(f.get("expected_outputs"));
			canon.append("subnet\t").append(f.get("order")).append('\t').append(f.get("id")).append('\t').append(f.get("name")).append('\t').append(f.get("type")).append('\t')
				.append(f.get("num_obs")).append('\t').append(f.get("expected_inputs")).append('\t').append(f.get("expected_outputs")).append('\t').append(f.get("seed")).append('\t')
				.append(f.get("dims")).append('\t').append(f.get("family_ranks")).append('\t')
				.append(f.get("loose_net_name")).append('\t').append(f.get("loose_subset_name")).append('\t').append(f.get("subset_sha256")).append('\t')
				.append(f.get("net_sha256")).append('\t').append(f.get("net_bytes"));
			if(expOut==4) canon.append('\t').append(f.get("output_names")).append('\t').append(f.get("output_units"))
				.append('\t').append(f.get("output_clamp_ratio"));
			canon.append('\n');
		}
		String newHash=sha256(canon.toString().getBytes(StandardCharsets.UTF_8));
		return text.replaceFirst("canonical_payload_sha256=[0-9a-f]+","canonical_payload_sha256="+newHash);
	}
	static String sha256(byte[] x) throws Exception { MessageDigest d=MessageDigest.getInstance("SHA-256"); byte[] y=d.digest(x); StringBuilder b=new StringBuilder(); for(byte v:y) b.append(String.format("%02x",v&255)); return b.toString(); }
	static void check(boolean ok,String what){ if(!ok) throw new AssertionError("FAIL: "+what); }
}
