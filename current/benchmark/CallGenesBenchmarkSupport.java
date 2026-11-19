package benchmark;

import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

import ddl.DDLFamilySearch.StructuralInputs;
import gff.GffLine;
import prok.ProkObject;
import prot.DigestSuffix;

/** Shared input, timing and grading contracts for the live caller benchmarks.
 * Historical diagnostic programs are not runtime dependencies of these helpers.
 * @author Sayu
 */
public final class CallGenesBenchmarkSupport{

	private CallGenesBenchmarkSupport(){}

	/** Whole-JVM CPU includes search workers and GC; current-thread CPU cannot replace it. */
	public static long processCpuNanos(){
		final java.lang.management.OperatingSystemMXBean bean=ProcessClock.BEAN;
		require(bean instanceof com.sun.management.OperatingSystemMXBean,
			"Process CPU clock is unavailable on this JVM");
		final long value=((com.sun.management.OperatingSystemMXBean)bean).getProcessCpuTime();
		require(value>=0, "Process CPU clock returned an unavailable value");
		return value;
	}

	/** Initialize the clock only when timing is requested, not when a grader reads its schema. */
	private static final class ProcessClock{
		private static final java.lang.management.OperatingSystemMXBean BEAN=ManagementFactory.getOperatingSystemMXBean();
	}

	/** Enable CDS grading and explicitly disable each inherited RNA grading flag. */
	public static void cdsOnly(){
		ProkObject.callCDS=true;
		ProkObject.calltRNA=ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=false;
	}

	/** Biological partial sides; absent and false are complete unless ranges state otherwise.
	 * Validate the partial flag first, including duplicate and malformed attributes.
	 */
	public static String partial(GffLine ref){
		ref.partial();
		final String p=CallGenesStopAudit.attribute(ref, "partial");
		final String low=CallGenesStopAudit.attribute(ref, "start_range");
		final String high=CallGenesStopAudit.attribute(ref, "end_range");
		final int l=range(low, ref.start), r=range(high, ref.stop);
		if(l<0 || r<0){return "unknown";}
		boolean left=l==1, right=r==1;
		if("true".equals(p)){
			if(!left && !right){return "unknown";}
		}else if(p!=null){
			final boolean pl=p.length()==2 && p.charAt(0)=='1', pr=p.length()==2 && p.charAt(1)=='1';
			if(low!=null && left!=pl || high!=null && right!=pr){return "unknown";}
			left=pl;
			right=pr;
		}
		if(ref.strand!=0 && ref.strand!=1){return "unknown";}
		final boolean five=ref.strand==0 ? left : right, three=ref.strand==0 ? right : left;
		return five && three ? "both" : five ? "five_prime" : three ? "three_prime" : "complete";
	}

	/** A singleton bound is complete; one open bound at the coordinate is partial. */
	private static int range(String text, int coordinate){
		if(text==null){return 0;}
		final int comma=text.indexOf(',');
		require(comma>0 && comma<text.length()-1 && text.indexOf(',', comma+1)<0,
			"Malformed partial range: "+text);
		final String a=text.substring(0, comma), b=text.substring(comma+1);
		final long low=a.equals(".") ? -1 : Long.parseLong(a), high=b.equals(".") ? -1 : Long.parseLong(b);
		require((a.equals(".") || low>0) && (b.equals(".") || high>0), "Invalid partial range: "+text);
		if(low==-1 && high==coordinate || high==-1 && low==coordinate){return 1;}
		return low==coordinate && high==coordinate ? 0 : -1;
	}

	/** Saved-reducer input policy. Structural mode establishes observed file stability,
	 * not byte identity: inherited pins are labels, never newly measured. Every declared
	 * input is checked before reading and again before publication.
	 */
	public static final class InputVerification{
		public InputVerification(String mode){
			require(mode==null || mode.equals("sha80") || mode.equals("structural"),
				"verification must be sha80 or structural");
			structural="structural".equals(mode);
		}

		public String mode(){return structural ? "structural" : "sha80";}

		public void pin(Path file, String label)throws Exception{
			require(label!=null && (label.matches("[0-9a-f]{20}") || structural && label.equals("unmeasured")),
				"Input requires sha80, or unmeasured with verification=structural");
			require(Files.isRegularFile(file, LinkOption.NOFOLLOW_LINKS) && !Files.isSymbolicLink(file),
				"Regular nonsymlink input required: "+file);
			final Path path=file.toAbsolutePath().normalize();
			final String prior=pins.get(path);
			require(prior==null || prior.equals(label), "Conflicting input identities: "+file);
			if(structural){
				stability.track(path);
				final Path real=path.toRealPath(), before=realPaths.putIfAbsent(path, real);
				require(before==null || before.equals(real), "Structural input changed: "+path);
			}else if(prior==null){
				require(label.equals(DigestSuffix.file(path.toString())), "Input pin mismatch: "+path);
			}
			pins.put(path, label);
		}

		public void verifyComplete()throws Exception{
			if(structural){
				stability.verify();
				for(Map.Entry<Path, Path> entry:realPaths.entrySet()){
					require(entry.getValue().equals(entry.getKey().toRealPath()), "Structural input changed: "+entry.getKey());
				}
			}else{
				for(Map.Entry<Path, String> entry:pins.entrySet()){
					require(entry.getValue().equals(DigestSuffix.file(entry.getKey().toString())), "Input changed: "+entry.getKey());
				}
			}
		}

		public final boolean structural;
		private final LinkedHashMap<Path, String> pins=new LinkedHashMap<Path, String>();
		private final StructuralInputs stability=new StructuralInputs();
		private final LinkedHashMap<Path, Path> realPaths=new LinkedHashMap<Path, Path>();
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Accepted saved-input column order, including retained historical identity labels. */
	public static final String INPUT_HEADER="genome\tfasta\tfasta_sha80\treference\treference_sha80\tnative\tnative_sha80"
		+"\tbaseline\tbaseline_sha80\tannotations\tannotations_sha80\tgrades\tgrades_sha80";

	/** Saved paired-caller call inventory; paths are data, never an instruction to rerun calling. */
	public static final String RETAINED_GFF_HEADER="genome\tnative_gff\tformula_gff";

	/** The native grader's eight reference- and query-relative counters, in emitted order. */
	public static final String[] GRADE_NAMES={"tp_start_reference", "tp_stop_reference", "fn_start", "fn_stop",
		"tp_start_query", "tp_stop_query", "fp_start", "fp_stop"};

	/** Public saved raw-score ABI. Historical private scorecards retain their frozen literals. */
	public static final String RAW_SCORE_HEADER="policy\tgenome\tannotation_or_called_length_band\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors\tunique_matched_keys\tduplicate_matching_occurrences\t"+String.join("\t",GRADE_NAMES)+"\n";
	/** Path-only saved-score input ABI; no sequence payload is read to emit this table. */
	public static final String SAVED_SCORE_INPUT_HEADER="policy\tgenome\tfasta\treference\tcalls\tgrades\tgrade_arm";
	/** Frozen raw controls carried into saved-output grading. */
	public static final String SAVED_SCORE_EXPECTED_HEADER="policy\tgenomes\treference_keys\tcorrect\tmissed\textras\tstop_errors\tstart_errors";
}
