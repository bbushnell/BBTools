package assemble;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import parse.Parse;
import shared.Shared;
import shared.Tools;

/** Optional final Hollow scaffolding for the public Tadpole single/multi-K entry points.
 * Original contigs and evidence reports remain in a unique sibling work directory.
 * @author Brian Bushnell, Fischl */
final class TadpoleScaffold {

	/** Removes scaffold options before native parsing; ordinary invocations take the unchanged path. */
	static boolean dispatch(final String[] input, final boolean forceMulti){
		assert(input!=null) : "Scaffolding dispatch follows native config expansion.";
		boolean seen=false, enabled=false;
		final ArrayList<String> core=new ArrayList<String>(), overrides=new ArrayList<String>();
		for(String arg : input){
			final String key=key(arg), value=value(arg);
			if(key.equals("scaffold")){seen=true; enabled=Parse.parseBoolean(value);}
			else if(key.startsWith("scaffold_")){
				seen=true;
				String option=key.substring(9);
				if(option.equals("threads")){option="t";}
				if(!Arrays.asList("k", "span", "tipbases", "tipmin", "tipmax", "slide", "minsupport",
						"gapratio", "flankcheck", "flankbases", "flankmin", "fillgaps", "t").contains(option)){
					throw new IllegalArgumentException("Tadpole owns scaffold paths; unsupported override: "+key);
				}
				overrides.add(option+"="+(value==null ? "t" : value));
			}else{core.add(arg);}
		}
		if(!seen){return false;}
		if(!enabled){TadpoleAssemblyPolicy.run(core.toArray(new String[0]), forceMulti); return true;}
		run(core, overrides, forceMulti);
		return true;
	}

	/** Builds the ordinary assembly first, then publishes only a successfully scaffolded result. */
	private static void run(final ArrayList<String> core, final ArrayList<String> overrides, final boolean forceMulti){
		assert(core!=null && overrides!=null) : "Dispatch separates native assembly arguments from postpass overrides.";
		final String[] original=core.toArray(new String[0]);
		if(!TadpoleAssemblyPolicy.isAssembly(original) || TadpoleGraph.requested(original)){
			throw new IllegalArgumentException("scaffold=t applies only to ordinary contig assembly, not correction, extension or graph-only mode.");
		}
		String in=null, mate=null, out=null;
		boolean overwrite=Tadpole.overwrite;
		final ArrayList<String> extras=new ArrayList<String>(), secondary=new ArrayList<String>();
		for(String arg : core){
			final String key=key(arg), value=value(arg);
			if(key.equals("in") || key.equals("in1")){in=value;}
			else if(key.equals("in2")){mate=value;}
			else if(key.equals("out") || key.equals("out1")){out=value;}
			else if(key.equals("extra") && present(value)){extras.addAll(Arrays.asList(value.split(",", -1)));}
			else if(key.equals("overwrite") || key.equals("ow")){overwrite=Parse.parseBoolean(value);}
			else if((key.equals("append") || key.equals("app")) && Parse.parseBoolean(value)){
				throw new IllegalArgumentException("Scaffolding cannot append to an existing assembly.");
			}else if(key.equals("out2") && present(value)){
				throw new IllegalArgumentException("Scaffolding publishes one FASTA; out2 is not supported.");
			}else if(present(value) && (key.startsWith("out") || key.equals("gfa") || key.equals("dot")
					|| key.equals("dump") || key.equals("hist") || key.equals("khist") || key.equals("ihist"))){
				secondary.add(value);
			}
		}
		if(!present(in) || !present(out) || out.indexOf(',')>=0 || out.indexOf('#')>=0
				|| FileFormat.isStdio(out) || !FileFormat.isFastaFile(out)){
			throw new IllegalArgumentException("scaffold=t requires named reads and one named out=FASTA (no output # placeholder).");
		}
		final ArrayList<String> reads=new ArrayList<String>();
		for(String path : in.split(",", -1)){
			if(path.indexOf('#')>=0 && !new File(path).exists()){
				final int pos=path.lastIndexOf('#');
				reads.add(path.substring(0, pos)+"1"+path.substring(pos+1));
				reads.add(path.substring(0, pos)+"2"+path.substring(pos+1));
			}else{reads.add(path);}
		}
		if(present(mate)){reads.addAll(Arrays.asList(mate.split(",", -1)));}
		reads.addAll(extras);
		// Native comma-delimited mate lists allow null placeholders for unpaired libraries.
		for(int i=reads.size()-1; i>=0; i--){if(reads.get(i).equalsIgnoreCase("null")){reads.remove(i);}}
		if(reads.isEmpty()){throw new IllegalArgumentException("Scaffolding requires at least one evidence file.");}
		for(String path : reads){
			if(path.isEmpty() || FileFormat.isStdio(path)){
				throw new IllegalArgumentException("Scaffolding rereads named evidence files; streams are not supported.");
			}
		}
		Tools.testInputFiles(false, true, reads);
		secondary.add(out);
		TadpoleGraph.checkPaths(reads.get(0), new ArrayList<String>(reads.subList(1, reads.size())), secondary.toArray(new String[0]));
		// Validate original publication against enabled model files before replacing out with a private intermediate.
		try(TadpoleAssemblyPolicy validation=new TadpoleAssemblyPolicy(original)){}
		final Path destination=Paths.get(out).toAbsolutePath().normalize();
		Path work=null;
		try{
			if(!Files.isDirectory(destination.getParent()) || Files.isDirectory(destination)
					|| (Files.exists(destination) && !overwrite)){
				throw new IOException("Output exists or parent is unavailable: "+destination);
			}
			work=Files.createTempDirectory(destination.getParent(), "tadpole-scaffold-");
			final Path contigs=work.resolve("contigs.fa"), scaffolds=work.resolve("scaffold-"+destination.getFileName());
			final ArrayList<String> assembly=new ArrayList<String>();
			for(String arg : core){if(!key(arg).equals("out") && !key(arg).equals("out1")){assembly.add(arg);}}
			assembly.add("out="+contigs);
			TadpoleAssemblyPolicy.run(assembly.toArray(new String[0]), forceMulti);
			if(!Files.isRegularFile(contigs)){throw new IOException("Assembly did not produce "+contigs);}
			if(Shared.AMINO_IN){throw new IllegalArgumentException("Hollow scaffolding requires nucleotide sequence.");}
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("k=224", "tipbases=192", "tipmin=10", "tipmax=160",
					"slide=0", "flankcheck=t", "flankbases=128", "flankmin=16", "minsupport=3", "gapratio=4",
					"fillgaps=f", "t="+Shared.threads()));
			args.addAll(overrides);
			args.addAll(Arrays.asList("join=t", "interleaved=f", "ow=f", "contigs="+contigs, "in="+reads.get(0),
					"out="+work.resolve("candidates.tsv"), "outjoins="+work.resolve("joins.tsv"), "outassembly="+scaffolds));
			if(reads.size()>1){args.add("extra="+String.join(",", reads.subList(1, reads.size())));}
			HollowBridge.main(args.toArray(new String[0]));
			if(overwrite){Files.move(scaffolds, destination, StandardCopyOption.REPLACE_EXISTING);}
			else{Files.move(scaffolds, destination);}
			System.err.println("Unscaffolded Tadpole contigs and link reports retained: "+work);
		}catch(IOException e){throw new RuntimeException("Scaffolding failed; intermediate evidence retained in "+work, e);}
	}

	/** Retains filename case and embedded equals signs. */
	private static String value(final String arg){final int p=arg.indexOf('='); return p<0 ? null : arg.substring(p+1);}
	/** Mirrors native option-name normalization without modifying path values. */
	private static String key(final String arg){
		final int p=arg.indexOf('=');
		String key=(p<0 ? arg : arg.substring(0, p)).toLowerCase(Locale.ROOT);
		while(key.startsWith("-")){key=key.substring(1);}
		return key;
	}
	private static boolean present(final String value){return value!=null && !value.isEmpty() && !value.equalsIgnoreCase("null");}
}
