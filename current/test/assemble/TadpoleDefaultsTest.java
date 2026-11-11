package test.assemble;

import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Compares direct construction and factory dispatch under an explicit defaults policy.
 * Each invocation constructs one short- or long-K assembler after setting known shared
 * values. No reads are loaded, workers started or output files opened by the assembler.
 * Run each case in a fresh JVM: the test deliberately changes process-wide settings.
 * This checks constructor configuration, not a complete BBMerge or Tadpole pipeline.
 * @author Brian Bushnell, Shinobu
 */
public final class TadpoleDefaultsTest{

	/** Static test entry point only. */
	private TadpoleDefaultsTest(){}

	/** Accepts in=existing.fa, k=31|63, factory=t|f and defaults=t|f. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TadpoleDefaultsTest.class, false);
		try{
			String input=null;
			int k=0;
			Boolean factory=null, defaults=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected named fixture options: "+arg);}
				if(pair[0].equalsIgnoreCase("in")){input=pair[1];}
				else if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(pair[1]);}
				else if(pair[0].equalsIgnoreCase("factory")){factory=Parse.parseBoolean(pair[1]);}
				else if(pair[0].equalsIgnoreCase("defaults")){defaults=Parse.parseBoolean(pair[1]);}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			if(input==null || input.length()==0 || (k!=31 && k!=63) || factory==null || defaults==null){
				throw new IllegalArgumentException("Required in=existing.fa k=31|63 factory=t|f defaults=t|f.");
			}
			run(input, k, factory.booleanValue(), defaults.booleanValue());
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Observes one constructor and checks the caller's requested policy independently of dispatch. */
	private static void run(final String input, final int k, final boolean factory, final boolean defaults){
		assert(k==31 || k==63) : "Fixture lengths select the two native table implementations.";
		seedSettings(k);
		final String[] args={"in="+input, "out=null", "k="+k, "mode=extend", "t=4",
				"prefilter=0", "prealloc=f", "ways=1", "initialsize=128", "ecc=f", "ecco=f",
				"merge=f", "markerrors=f", "hpindel=f", "fixindels=f", "hollow=f",
				"hollowbranchstop=f", "hollowbranchrollback=f", "hollowrescue=f",
				"el=0", "er=1", "showstats=f"};
		final Tadpole tad=factory ? Tadpole.makeTadpole(args, defaults) :
				k==31 ? new Tadpole1(args, defaults) : new Tadpole2(args, defaults);
		try{
			check(tad.getClass()==(k==31 ? Tadpole1.class : Tadpole2.class), "Wrong table implementation.");
			final AbstractKmerTableSet table=tad.tables();
			final float minProb=k==31 ? ((KmerTableSet)table).minProb : ((KmerTableSetU)table).minProb;
			final String observed="k="+k+" factory="+factory+" defaults="+defaults+
					" implementation="+tad.getClass().getSimpleName()+" threads="+Shared.threads()+
					" zip="+ReadWrite.ZIPLEVEL+" pigz="+ReadWrite.USE_PIGZ+" unpigz="+ReadWrite.USE_UNPIGZ+
					" bf1="+ByteFile.FORCE_MODE_BF1+" bf2="+ByteFile.FORCE_MODE_BF2+
					" bf3="+ByteFile.FORCE_MODE_BF3+" bf4="+ByteFile.FORCE_MODE_BF4+
					" global_minprob="+AbstractKmerTableSet.defaultMinprob+" table_minprob="+minProb;
			System.out.println("DEFAULTS_OBSERVATION "+observed);
			final double expectedMinProb=defaults ? 0.5 : 0.125;
			check(Shared.threads()==4 && !ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF3 &&
					!ByteFile.FORCE_MODE_BF4, "Unexpected fixture thread/reader setup: "+observed);
			check(ReadWrite.ZIPLEVEL==(defaults ? 2 : 9) && ReadWrite.USE_PIGZ==defaults &&
					ReadWrite.USE_UNPIGZ==defaults && ByteFile.FORCE_MODE_BF2==defaults &&
					AbstractKmerTableSet.defaultMinprob==expectedMinProb && minProb==(float)expectedMinProb,
					"DEFAULTS_CONTRACT_MISMATCH "+observed);
			System.out.println("DEFAULTS_CONTRACT_PASS");
		}finally{tad.clearData();}
	}

	/** Seeds values distinct from Tadpole defaults and aligns representation state for both entry paths. */
	private static void seedSettings(final int k){
		Kmer.PACKED=false;
		AbstractKmerTableSet.MASK_CORE=k<=31;
		Shared.COMMAND_LINE=null;
		Shared.setThreads(4);
		ReadWrite.ZIPLEVEL=9;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		ByteFile.FORCE_MODE_BF1=ByteFile.FORCE_MODE_BF2=false;
		ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
		AbstractKmerTableSet.defaultMinprob=0.125;
		assert(Shared.threads()==4) : "Four configured threads exercise the constructor's BF2 default branch.";
	}

	/** Contract checks remain active independently of the JVM assertion setting. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
