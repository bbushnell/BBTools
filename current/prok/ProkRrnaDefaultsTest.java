package prok;

import java.util.ArrayList;

/** Exercises the actual parser/registration path in separate JVMs.
 * @author Brian Bushnell, Raiden */
public final class ProkRrnaDefaultsTest {
	public static void main(String[] args){
		if(args.length!=2){throw new IllegalArgumentException("Expected test mode and input FASTA");}
		final String mode=args[0];final ArrayList<String> flags=new ArrayList<String>();
		flags.add("in="+args[1]);flags.add("t=1");flags.add("taxonomy=f");
		flags.add("trnaalign=f");flags.add("trnaboundarynet=f");
		int expected=ProkRrnaDefaults.S16|ProkRrnaDefaults.S23|ProkRrnaDefaults.S5;
		if(mode.equals("pgm")){flags.add("rrna17=f");expected=0;}
		else if(mode.equals("disable16")){flags.add("call16s=f");expected=ProkRrnaDefaults.S23|ProkRrnaDefaults.S5;}
		else if(mode.equals("disable23")){flags.add("call23s=f");expected=ProkRrnaDefaults.S16|ProkRrnaDefaults.S5;}
		else if(mode.equals("disable5")){flags.add("call5s=f");expected=ProkRrnaDefaults.S16|ProkRrnaDefaults.S23;}
		else if(mode.equals("endpoints-off")){flags.add("rrna5sendpointnn=f");}
		else if(mode.equals("disableall")){flags.add("call16s=f");flags.add("call23s=f");flags.add("call5s=f");expected=0;}
		else if(mode.equals("explicit23")){flags.add("rrna17=23s");flags.add("call23s=f");expected=ProkRrnaDefaults.S23;}
		else if(mode.equals("default-controls")){
			flags.add("rrna16sseedminhits=2");flags.add("rrna23sidpass=0.60");
			flags.add("rrna16swindow=containment");flags.add("rrna23swindow=legacy");
			flags.add("rrna16spartialminlen=1200");flags.add("rrna23spartialminlen=1800");
		}
		else if(!mode.equals("default")){throw new IllegalArgumentException("Unknown test mode "+mode);}
		new CallGenes(flags.toArray(new String[0]));
		int observed=0;
		for(NcrnaFamily family:GeneCaller.ncrnaFamilies){
			if(family.name.equals("16S")){observed|=ProkRrnaDefaults.S16;}
			else if(family.name.equals("23S")){observed|=ProkRrnaDefaults.S23;}
			else if(family.name.equals("prok5S")){
				observed|=ProkRrnaDefaults.S5;
				assert((family.rrnaEndpointFeatures!=null)==!mode.equals("endpoints-off"))
					: "Seeded5S must load its bundled endpoint nets unless explicitly disabled";
			}
			else{throw new AssertionError("No unrelated RNA family may be enabled by prokaryotic defaults: "+family.name);}
		}
		assert(observed==expected) : "Actual parser/loader family set differs for "+mode+": "+observed+" != "+expected;
		assert(ProkObject.call16S==(mode.equals("pgm") || mode.equals("explicit23"))) : "PGM16S must be restored by rrna17=f or preserved by a23S-only profile";
		assert(ProkObject.call23S==mode.equals("pgm")) : "PGM23S must be restored only when enabled and not replaced by the chosen seeded method";
		assert((observed&ProkRrnaDefaults.S16)==0 || !ProkObject.call16S) : "A seeded16S default must disable duplicate PGM16S";
		assert((observed&ProkRrnaDefaults.S23)==0 || !ProkObject.call23S) : "A seeded23S default must disable duplicate PGM23S";
		assert(ProkObject.call5S==(mode.equals("pgm") || mode.equals("explicit23")))
			: "PGM5S returns with explicit PGM or an unrelated explicit profile; seeded defaults must not duplicate it";
		assert(ProkObject.callCDS && ProkObject.calltRNA && !ProkObject.call18S) : "Default method selection must preserve CDS/tRNA and the disabled18S gate";
		assert(ProkRrnaDefaults.mask(ProkRrnaDefaults.Method.PGM,ProkRrnaDefaults.Method.PGM,ProkRrnaDefaults.Method.RRNA17,true,true,true)==ProkRrnaDefaults.S5)
			: "A future5S-only constant change must not also enable16S/23S";
		assert(ProkRrnaDefaults.mask(ProkRrnaDefaults.Method.RRNA17,ProkRrnaDefaults.Method.RRNA17,ProkRrnaDefaults.Method.RRNA17,false,false,false)==0)
			: "Explicit family disabling must dominate every default method";
		System.out.println("PROK_RRNA_DEFAULTS_PASS mode="+mode+" mask="+observed);
	}
}
