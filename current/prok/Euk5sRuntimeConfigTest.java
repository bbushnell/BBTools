package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

/** Fresh-JVM CLI fixtures; build producer supplies a test-only seven-model/A bundle.
 * @author Raiden
 */
public final class Euk5sRuntimeConfigTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1,"Supply one fresh-JVM case: on, off, off_endpoint, off_cutoffs, off_enabled, off_path, scalar, ablation, missing, table");
		final String mode=args[0];final Path dir=Files.createTempDirectory("euk5s_runtime_");
		check(mode.matches("on|off|off_endpoint|off_cutoffs|off_enabled|off_path|scalar|ablation|missing|table"),"Unknown fixture case: "+mode);
		final Path input=dir.resolve("input.fa");Files.write(input,">fixture\nACGTACGTACGTACGTACGT\n".getBytes(StandardCharsets.US_ASCII));
		final ArrayList<String> flags=new ArrayList<String>();
		for(String flag:new String[]{"in="+input,"out=null","taxonomy=f","trnaboundarynet=f","t=1"}){flags.add(flag);}
		final boolean off=mode.startsWith("off");flags.add("euk5s="+(off?"f":"t"));
		final boolean failure=mode.equals("missing") || mode.equals("off_enabled") || mode.equals("off_path");
		if(mode.equals("scalar")){flags.add("ncrnafamily=euk5S");flags.add("ncrnaidpass=.71");flags.add("ncrnaidborderline=.70");}
		else if(mode.equals("ablation")){flags.add("euk5sendpoint=f");}
		else if(mode.equals("missing")){flags.add("euk5sendpointnets="+dir.resolve("missing"));}
		else if(mode.equals("off_endpoint")){flags.add("euk5sendpoint=f");}
		else if(mode.equals("off_cutoffs")){flags.add("euk5smodelcutoffs=f");}
		else if(mode.equals("off_enabled")){flags.add("euk5sendpoint=t");}
		else if(mode.equals("off_path")){flags.add("euk5sendpointnets="+dir.resolve("missing"));}
		else if(mode.equals("table")){
			final String[] names={"euk5S_universal","euk5S_fungi","euk5S_plant","euk5S_animal","euk5S_protist","euk5S_residual_protozoa","euk5S_residual_invertebrate"};
			final structures.ByteBuilder b=new structures.ByteBuilder("family\tmodel\tidpass\tidborderline\n");
			for(int i=names.length-1;i>=0;i--){b.append("euk5S\t").append(names[i]).append("\t.72\t.71\n");}
			final Path table=dir.resolve("cutoffs.tsv");Files.write(table,b.toBytes());
			flags.add("euk5smodelcutoffs="+table);flags.add("ncrnafamily=euk5S");flags.add("ncrnaidpass=.75");flags.add("ncrnaidborderline=.74");
		}else{check(off || mode.equals("on"),"Unknown fixture case: "+mode);}
		if(off){
			// All euk5S resource resolution must be skipped for an empty family list.
			final Euk5sRuntimeConfig config=new Euk5sRuntimeConfig();config.apply(new ArrayList<NcrnaFamily>());
		}
		try{new CallGenes(flags.toArray(new String[0]));}
		catch(RuntimeException e){if(failure){
			check(e.toString().contains("missing") || e.toString().contains("MANIFEST") || e.toString().contains("runtime options require"),"Resource failure must name its missing dependency: "+e);
			System.out.println("PASS Euk5sRuntimeConfigTest "+mode);return;
		}throw e;}
		check(!failure,"Missing endpoint resources or disabled-family dependencies must fail before calling");
		if(off){check(GeneCaller.ncrnaFamilies.isEmpty(),"euk5s=f must not register any new family");}
		else{
			check(GeneCaller.ncrnaFamilies.size()==1,"Ordinary CLI must independently register exactly euk5S");
			final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(0);
			check(family.library.length==7 && family.indexTopN==7 && family.seedMinHits==2 && family.kLong==17,"Ordinary CLI must install the seven-model/two-hit configuration");
			if(mode.equals("scalar")){check(family.modelThresholds==null && family.idPass==.71f && family.idBorderline==.70f,"Scalar sweep must reach the caller instead of being masked by the default table");}
			else{
				check(family.modelThresholds!=null,"Ordinary CLI must load a model-bound cutoff table");
				for(int i=0;i<7;i++){check(family.modelThresholds.pass(i)==(mode.equals("table")?.72f:.68f)
					&& family.modelThresholds.borderline(i)==(mode.equals("table")?.71f:.68f),"Explicit tables bind by model and override scalar controls");}
			}
			if(mode.equals("ablation")){check(family.rrnaEndpointFeatures==null,"Endpoint ablation must retain a null inference hook");}
			else{check(family.rrnaEndpointFeatures!=null && family.rrnaEndpointFeatures.fiveNets.length==7
				&& family.rrnaEndpointFeatures.threeNets.length==7,"Ordinary CLI must install both endpoint networks for every model");}
			final NcrnaScavenger caller=RrnaEndpointCallerFeaturesTest.scavenger(family);
			check(caller!=null,"Runtime resources must also pass the actual worker configuration path");
		}
		System.out.println("PASS Euk5sRuntimeConfigTest "+mode);
	}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
