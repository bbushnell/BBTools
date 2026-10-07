package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import cardinality.DynamicDemiLog;
import clade.Clade;
import clade.SendClade;
import stream.Read;
import tracker.EntropyTracker;

/**
 * Exercises the server-response boundary independently of network availability and reference assets.
 * A malformed response must never become a valid unknown feature row.
 * @author Yoimiya
 */
public final class MagQCAssemblyInputTest {

	/** Runs classified, partial, no-hit, identity and malformed-envelope cases. */
	public static void main(String[] args) throws Exception{
		check(args.length==0, "This fixture accepts no arguments");
		checkOverrides();
		checkAmbiguousTransport();
		checkSketchSession();
		checkNormalAcknowledgment();
		checkTransportFailure();
		check(parse("#Query1\n").status.equals("unknown"), "Valid no-hit envelope must remain unknown");
		final MagQCAssemblyInput.Taxonomy classified=parse(row("d__Bacteria;p__Bacillota;g__Bacillus", "d:0.0;p:0.0"));
		check(classified.status.equals("classified") && classified.domain.equals("Bacteria") &&
			classified.phylum.equals("Bacillota"), "C1 uses lineage, not per-rank confidence, for features");
		check(parse(row("d__Archaea", ".")).status.equals("partial"), "Domain-only lineage must be partial");
		check(parse(row("NA", ".")).status.equals("unknown"), "Valid result with no lineage must be unknown");
		check(parse(row("d:Bacteria;p:Bacillota", "d:0.0;p:0.0")).phylum.equals("Bacillota"),
			"CallGenes also accepts colon-prefixed lineage without confusing it with confidence");
		check(parse(row("d__Bacteria;sk__Bacteria;p__Bacillota;ss__a;ss__b;st__c", ".")).phylum.equals("Bacillota"),
			"Repeated subspecies and multi-letter prefixes emitted by Clade.lineage are valid");
		// Optional SSU and sketch metrics precede lineage; reference cardinality is the sixth DDL field.
		check(parse(row("d__Bacteria;p__Bacillota", ".").replace("\t0.9\t0.8\t0.7", "\t0.8\t0.7")).status.equals("classified"),
			"Lineage extraction must not shift onto reference cardinality when SSU is absent");
		check(parse(row("d__Bacteria;p__Bacillota", ".").replace("\t0.9\t0.8\t0.7\t0.6\t0.5\t100\t200", "")).status.equals("classified"),
			"Machine results without either optional SSU or DDL columns must retain lineage");
		reject(null); reject(""); reject("\n"); reject("Internal server error\n");
		reject("No valid clades found in request\n"); reject("#Query1\nInternal server error\n");
		reject("#Query2\n"); reject("#Query1\n#Query2\n");
		final String good=row("d__Bacteria;p__Bacillota", ".");
		reject(good+good); reject(good+good.substring(good.indexOf('\n')+1));
		reject(good.replace("magqc_bin", "wrong_query"));
		reject(good.replace("\t1000\t2\t", "\t999\t2\t"));
		reject(good.replace("\t1000\t2\t", "\t1000\t3\t"));
		reject(good.replace("0.500", "NaN"));
		reject(good.replace("d__Bacteria;p__Bacillota", "garbled taxonomy"));
		reject(row("p__Bacillota", "."));
		reject(row("d__Bacteria;d__Archaea;p__Bacillota", "."));
		reject(row("d__Bacteria;p__Bacillota;p__Other", "."));
		System.out.println("MagQCAssemblyInputTest PASS: explicit taxonomy, C1 lineage, no-hit, envelope and query identity");
	}

	/** An older server must not silently satisfy a request for normal candidate search. */
	private static void checkNormalAcknowledgment(){
		final String body=row("d__Bacteria;p__Bacillota", ".");
		final String marker="#QuickCladeNormalSearch";
		for(String ending:new String[]{"\n", "\r\n"}){
			final String stripped=SendClade.requireNormalAck(marker+ending+body);
			check(body.equals(stripped) && parse(stripped).phylum.equals("Bacillota"),
				"Normal acknowledgment must preserve the complete machine response for taxonomy parsing");
		}
		for(String invalid:new String[]{null, "", body, marker, "#QuickCladeNormalSearchExtra\n"+body, body+marker+"\n"}){
			boolean rejected=false;
			try{SendClade.requireNormalAck(invalid);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Absent, malformed or misplaced normal acknowledgment must fail before taxonomy is accepted");
		}
	}

	/** A malformed endpoint must release the throttle slot before the next batch request. */
	private static void checkTransportFailure() throws Exception{
		final java.lang.reflect.Field field=SendClade.class.getDeclaredField("concurrency");
		field.setAccessible(true);
		final java.util.concurrent.atomic.AtomicInteger slots=(java.util.concurrent.atomic.AtomicInteger)field.get(null);
		final int before=slots.get();
		check(!SendClade.sync && before==0, "Transport fixture requires an idle asynchronous throttle");
		// Invalid URI syntax fails before any network access in ServerTools.
		final String result=SendClade.sendMessage(new byte[]{1}, "not a URI", false);
		check(result==null, "Malformed endpoint became a valid classifier reply");
		check(slots.get()==before, "Malformed endpoint leaked a concurrency slot and can deadlock later requests");
		System.out.println("Transport failure PASS: malformed endpoint releases throttle slot without network access");
	}

	/** Compares actual serialized sketches and tests lifecycle exclusion with overlapping workers. */
	private static void checkSketchSession() throws Exception{
		final ArrayList<Read> contigs=new ArrayList<Read>();
		final byte[] bases=new byte[20000];
		final byte[] alphabet="ACGTN".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final java.util.Random random=new java.util.Random(61);
		for(int i=0; i<bases.length; i++){bases[i]=alphabet[random.nextInt(alphabet.length)];}
		contigs.add(new Read(bases, null, "arbitrary_header_tid_999", 0));
		final boolean oldDdl=Clade.MAKE_DDLS;
		final int oldK=Clade.DDL_K, oldBuckets=Clade.DDL_BUCKETS, oldExponent=DynamicDemiLog.exponentBits();
		final long oldSeed=Clade.DDL_SEED;
		final byte[] expected;
		try{
			// Independent copy of the pre-session C1 query construction, without HTTP.
			Clade.MAKE_DDLS=true; Clade.DDL_K=25; Clade.DDL_BUCKETS=32768; Clade.DDL_SEED=12345L;
			DynamicDemiLog.setExponent(5); bin.AdjustEntropy.load(4, 150);
			final Clade query=new Clade(0, 0, MagQCAssemblyInput.QUERY_NAME);
			final EntropyTracker entropy=new EntropyTracker(4, 150, false);
			for(Read read:contigs){query.add(read.bases, entropy);}
			query.finish();
			final ArrayList<Clade> queries=new ArrayList<Clade>(); queries.add(query);
			expected=SendClade.toMessage(queries, true, 1, false, false, 1, 1);
		}finally{
			Clade.MAKE_DDLS=oldDdl; Clade.DDL_K=oldK; Clade.DDL_BUCKETS=oldBuckets; Clade.DDL_SEED=oldSeed;
			DynamicDemiLog.setExponent(oldExponent);
		}
		final MagQCAssemblyInput.SketchSession session=MagQCAssemblyInput.openSketchSession();
		try{
			check(Arrays.equals(expected, session.request(contigs).bytes), "Session changes the legacy C1 request bytes");
			boolean rejected=false;
			try{MagQCAssemblyInput.openSketchSession();}catch(IllegalStateException e){rejected=true;}
			check(rejected, "Nested sessions must not restore another worker's sketch globals");
			rejected=false;
			try{MagQCAssemblyInput.classify(contigs, "refseq");}catch(IllegalStateException e){rejected=true;}
			check(rejected, "Legacy classifier must reject setup changes before any HTTP call");
			final CountDownLatch entered=new CountDownLatch(2), release=new CountDownLatch(1);
			final ArrayList<Read> held=new HeldContigs(contigs, entered, release);
			final Throwable[] errors=new Throwable[2];
			final Thread[] threads=new Thread[2];
			for(int i=0; i<threads.length; i++){
				final int index=i;
				threads[i]=new Thread(new Runnable(){
					@Override public void run(){
						try{
							final MagQCAssemblyInput.SketchRequest request=session.request(held);
							check(Arrays.equals(expected, request.bytes), "Parallel C1 wire bytes differ");
							check(request.bases==20000 && request.contigs==1, "Echoed counts omit undefined bases");
						}catch(Throwable t){errors[index]=t;}
					}
				}, "taxonomy-fixture-"+i);
				threads[i].start();
			}
			try{
				check(entered.await(10, TimeUnit.SECONDS), "Whole-query locking prevents concurrent sketch entry");
				rejected=false;
				try{session.close();}catch(IllegalStateException e){rejected=true;}
				check(rejected, "Closing active workers must not restore their globals");
			}finally{
				release.countDown();
				for(Thread thread:threads){thread.join();}
			}
			for(Throwable error:errors){if(error!=null){throw new AssertionError("Parallel sketch failed", error);}}
			check(session.peakWorkers()==2, "Fixture must exercise two overlapping sketch workers");
			rejected=false;
			try{session.request(new ArrayList<Read>());}catch(IllegalArgumentException e){rejected=true;}
			check(rejected && Arrays.equals(expected, session.request(contigs).bytes), "Bad input must release session accounting");
		}finally{session.close();}
		check(Clade.MAKE_DDLS==oldDdl && Clade.DDL_K==oldK && Clade.DDL_BUCKETS==oldBuckets &&
			Clade.DDL_SEED==oldSeed && DynamicDemiLog.exponentBits()==oldExponent, "Prior sketch settings were not restored");
		boolean rejected=false;
		try{session.request(contigs);}catch(IllegalStateException e){rejected=true;}
		check(rejected, "Closed session must not build with arbitrary restored globals");
		try(MagQCAssemblyInput.SketchSession normal=MagQCAssemblyInput.openSketchSession(true)){
			final byte[] normalBytes=normal.request(contigs).bytes;
			final String legacyText=new String(expected, java.nio.charset.StandardCharsets.UTF_8);
			final String normalText=new String(normalBytes, java.nio.charset.StandardCharsets.UTF_8);
			check(normalText.indexOf('\n')>0 && legacyText.indexOf('\n')>0,
				"Both search policies must serialize a request header followed by C1 sketches");
			check(!normalText.substring(0, normalText.indexOf('\n')).equals(legacyText.substring(0, legacyText.indexOf('\n'))),
				"Normal search must request a distinct policy instead of silently sending legacy bytes");
			check(normalText.substring(normalText.indexOf('\n')).equals(legacyText.substring(legacyText.indexOf('\n'))),
				"Search policy must not change the C1 query sketch recipe or serialized query data");
		}
		System.out.println("Sketch session PASS: legacy wire parity, overlapping workers, lifecycle exclusion and failure recovery");
	}

	/** Pauses only the fixture's iteration so the owner can inspect a deterministic active overlap. */
	private static final class HeldContigs extends ArrayList<Read>{
		HeldContigs(ArrayList<Read> source, CountDownLatch entered_, CountDownLatch release_){
			super(source); entered=entered_; release=release_;
		}
		@Override public Iterator<Read> iterator(){
			entered.countDown();
			try{check(release.await(10, TimeUnit.SECONDS), "Fixture owner failed to release workers");}
			catch(InterruptedException e){Thread.currentThread().interrupt(); throw new RuntimeException(e);}
			return super.iterator();
		}
		private final CountDownLatch entered, release;
		private static final long serialVersionUID=1L;
	}

	/** Offline taxonomy is explicit, domain-bound and identified separately in provenance. */
	private static void checkOverrides(){
		final java.util.HashMap<String,String> options=new java.util.HashMap<String,String>();
		check(MagQCAssemblyInput.overrideTaxonomy(options)==null, "Default still requires QuickClade classification");
		options.put("taxphylum", "Bacillota");
		boolean rejected=false;
		try{MagQCAssemblyInput.overrideTaxonomy(options);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A phylum must never cause the domain feature to be guessed");
		options.put("taxdomain", "bacteria");
		final MagQCAssemblyInput.Taxonomy full=MagQCAssemblyInput.overrideTaxonomy(options);
		check(full.domain.equals("Bacteria") && full.phylum.equals("Bacillota") && full.status.equals("classified"),
			"Explicit taxonomy must retain the supplied phylum with canonical domain spelling");
		options.remove("taxphylum"); options.put("taxdomain", "Archaea");
		check(MagQCAssemblyInput.overrideTaxonomy(options).status.equals("partial"), "Domain-only override must remain partial");
		options.put("fasta", "fixture.fa");
		final structures.ByteBuilder provenance=new structures.ByteBuilder();
		MagQCAssemblyInput.appendProvenance(provenance, options);
		check(provenance.toString().contains("#taxonomy_source\tuser-supplied\n"), "Report must distinguish override from server prediction");
		for(String invalid:new String[]{"", "Bacillota\tgarbage", "Bacillota\nextra", " Bacillota"}){
			options.put("taxphylum", invalid); rejected=false;
			try{MagQCAssemblyInput.overrideTaxonomy(options);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Empty or malformed phylum override must fail before output");
		}
		options.remove("taxphylum"); options.put("taxdomain", "Eukaryota"); rejected=false;
		try{MagQCAssemblyInput.overrideTaxonomy(options);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Prokaryote override requires an explicit supported domain");
	}

	/** Native serialization proves server Q_Bases includes the fifth, undefined-base bucket. */
	private static void checkAmbiguousTransport(){
		final boolean old=clade.Clade.MAKE_DDLS;
		try{
			clade.Clade.MAKE_DDLS=false;
			bin.AdjustEntropy.load(4, 150);
			final clade.Clade query=new clade.Clade(0, 0, MagQCAssemblyInput.QUERY_NAME);
			query.add("ACGTNNACGT".getBytes(java.nio.charset.StandardCharsets.US_ASCII),
				new tracker.EntropyTracker(4, 150, false));
			query.finish();
			final java.util.ArrayList<byte[]> lines=new java.util.ArrayList<byte[]>();
			for(String line:query.toBytes(null).toString().split("\n")){
				lines.add(line.getBytes(java.nio.charset.StandardCharsets.UTF_8));
			}
			final clade.Clade restored=clade.Clade.parseCladeFlex(lines, new parse.LineParser1('\t'));
			check(query.bases==10 && restored.bases==10 && restored.bases==query.monomerSum(),
				"Clade transport must include Ns: query="+query.bases+" restored="+restored.bases+
				" monomers="+query.monomerSum());
			final String response=row("d__Bacteria;p__Bacillota", ".").replace("\t1000\t2\t", "\t10\t1\t");
			check(MagQCAssemblyInput.parseResponse(response, query.monomerSum(), 1).status.equals("classified"),
				"An N-containing assembly must accept its correctly echoed transport counts");
		}finally{clade.Clade.MAKE_DDLS=old;}
	}

	/** Literal default machine columns follow Comparison.appendResultMachine(false). */
	private static String row(String lineage, String confidence){
		return "#Query1\nmagqc_bin\t0.500\t1000\t2\treference\t123\t0.501\t2000\t3\tspecies"+
			"\t0.001\t0.01\t0.02\t0.03\t0.04\t0.05\t0.06"+
			"\t0.9\t0.8\t0.7\t0.6\t0.5\t100\t200\t"+lineage+"\tdomain\t"+confidence+"\n";
	}

	/** Parses against the submitted query's independently known counts. */
	private static MagQCAssemblyInput.Taxonomy parse(String response){
		return MagQCAssemblyInput.parseResponse(response, 1000, 2);
	}

	/** A failure is required even with JVM assertions disabled. */
	private static void reject(String response){
		try{parse(response);}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed classifier response was accepted: "+response);
	}

	/** States the consequence of each response-contract check. */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
