package test.jgi;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;

import jgi.KmerCountExact;
import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import shared.Shared;

/** Checks CLI thread selection and unchanged default/public-field behavior after construction.
 * Each case needs a fresh JVM with the specified ActiveProcessorCount and no scheduler
 * environment override of Shared.LOGICAL_PROCESSORS. No loading or table allocation runs.
 * @author Sucrose
 */
public final class KmerCountExactThreadsTest{

	/** Runs one literal CLI/default/API case under either assertion mode. */
	public static void main(String[] args) throws ReflectiveOperationException{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=<thread case> argument.");
		}
		final String mode=args[0].substring(5);
		final int startup, expectedShared, expectedTool;
		final String option;
		final boolean api;
		switch(mode){
			case "threads-one": startup=2; expectedShared=expectedTool=1; option="threads=1"; api=false; break;
			case "t-two": startup=1; expectedShared=expectedTool=2; option="t=2"; api=false; break;
			case "threads-two-control": startup=2; expectedShared=expectedTool=2; option="threads=2"; api=false; break;
			case "t-one-control": startup=1; expectedShared=expectedTool=1; option="t=1"; api=false; break;
			case "default-control": startup=1; expectedShared=expectedTool=1; option=null; api=false; break;
			case "api-control": startup=1; expectedShared=1; expectedTool=2; option=null; api=true; break;
			case "auto-control": startup=2; expectedShared=expectedTool=2; option="threads=auto"; api=false; break;
			default: throw new IllegalArgumentException("Unknown thread case: "+mode);
		}
		//Read the real class-load default before changing any shared or tool thread field.
		final int initialShared=Shared.threads(), initialTool=KmerCountExact.THREADS;
		check(Runtime.getRuntime().availableProcessors()==startup && Shared.LOGICAL_PROCESSORS==startup
			&& initialShared==startup && initialTool==startup, mode,
			"Fixture needs startup processors="+startup+", available="+Runtime.getRuntime().availableProcessors()+
			", logical="+Shared.LOGICAL_PROCESSORS+", shared="+initialShared+", tool="+initialTool);
		KmerCountExact.DISPLAY_PROGRESS=false;
		AbstractKmerTableSet.DISPLAY_PROGRESS=false;
		if(api){KmerCountExact.THREADS=2;}
		final ArrayList<String> options=new ArrayList<String>(Arrays.asList(
			"in=/dev/null", "k=3", "ways=1", "prefilter=0", "prealloc=f", "onepass=f"));
		if(option!=null){options.add(option);}
		final KmerCountExact tool=new KmerCountExact(options.toArray(new String[options.size()]));
		checkUnallocated(tool, mode);
		check(Shared.threads()==expectedShared, mode,
			"Shared parser must select requested threads: expected="+expectedShared+", actual="+Shared.threads());
		check(KmerCountExact.THREADS==expectedTool, mode,
			"Tool thread field must match the selected CLI/default/API contract: expected="+expectedTool+
			", actual="+KmerCountExact.THREADS+", shared="+Shared.threads());
		System.out.println("PASS kmer-count-threads case="+mode+" startup="+startup+
			" shared="+Shared.threads()+" tool="+KmerCountExact.THREADS);
	}

	/** Observes the existing private table-set reference without modifying production state. */
	private static void checkUnallocated(final KmerCountExact tool, final String mode) throws ReflectiveOperationException{
		final Field field=KmerCountExact.class.getDeclaredField("tables");
		field.setAccessible(true);
		final AbstractKmerTableSet tables=(AbstractKmerTableSet)field.get(tool);
		check(tables instanceof KmerTableSet && tables.kbig()==3 && tables.ways()==1
			&& tables.tables()==null && tables.prefilterArray==null && !tables.prefilter && !tables.onePass
			&& tables.in1.size()==1 && tables.in1.get(0).equals("/dev/null") && tables.in2.isEmpty()
			&& tables.extra.isEmpty() && tables.readsIn==0 && tables.basesIn==0 && tables.kmersIn==0
			&& !tables.errorState, mode,
			"Fixture must stop after construction with no loaded reads, table arrays or prefilter");
	}

	/** Keeps the parser and tool-field oracles active with Java assertions disabled. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR119 case="+mode+": "+message);}
	}
}
