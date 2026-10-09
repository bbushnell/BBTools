package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;

/** Backfills native tip entropy without rereading reads or rebuilding assemblies. @author Fischl */
final class FusionEntropyCorpus {

	/** Builds an immutable109-input derivative and group-disjoint native trainer tables. */
	static void run(final String[] args){
		final Parser parser=new Parser();
		String root=null;
		long seed=350195;
		int expected=1000;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase(java.util.Locale.ROOT);
			while(a.startsWith("-")){a=a.substring(1);}
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("root")){root=b;}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("expected")){expected=Integer.parseInt(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown entropy corpus option: "+arg);}
		}
		if(parser.in1==null || parser.out1==null || root==null || expected<1){
			throw new IllegalArgumentException("Entropy mode requires in=samples.tsv root=corpus_root out=new_directory.");
		}
		if(!Tools.testInputFiles(false, true, parser.in1)){throw new IllegalArgumentException("Missing sample manifest.");}
		final ArrayList<Sample> samples=readSamples(parser.in1, expected);
		assignSplits(samples, seed);
		// Preflight all source paths before any derivative output is opened.
		for(Sample sample : samples){
			final String dir=root+"/data/"+sample.id;
			if(!Tools.testInputFiles(false, true, dir+"/joins.vectors.tsv", dir+"/joins.trace.tsv", dir+"/labels.tsv")){
				throw new IllegalArgumentException("Incomplete source corpus: "+sample.id);
			}
		}
		final File out=new File(parser.out1);
		if(out.exists() || !out.mkdir()){throw new IllegalArgumentException("Output directory must be new: "+out);}
		final ByteStreamWriter[] writers=new ByteStreamWriter[4];
		final long[] totals=new long[5];
		boolean error=false;
		try{
			for(int i=0; i<3; i++){
				writers[i]=open(out+"/"+SPLITS[i]+".tsv");
				writers[i].print("#dims\t109\t1\n");
			}
			writers[3]=open(out+"/census.tsv");
			writers[3].print("genome\tgroup\tsplit\tcandidates\tsupported\tcontradicted\tunresolved\tunavailable_tips\n");
			final FusionTipEntropy entropy=new FusionTipEntropy();
			final ByteBuilder row=new ByteBuilder();
			for(Sample sample : samples){
				final long[] counts=convert(sample, root, out.toString(), entropy, writers[sample.split]);
				row.clear().append(sample.id).tab().append(sample.group).tab().append(SPLITS[sample.split]);
				for(int i=0; i<counts.length; i++){row.tab().append(counts[i]); totals[i]+=counts[i];}
				writers[3].print(row.nl());
			}
		}finally{
			for(ByteStreamWriter writer : writers){if(writer!=null){error=writer.poisonAndWait() | error;}}
		}
		if(error){throw new IllegalStateException("Entropy corpus output failed.");}
		System.err.println("FUSION_ENTROPY_CORPUS_PASS genomes="+samples.size()+" rows="+totals[0]+
				" supported="+totals[1]+" contradicted="+totals[2]+" unresolved="+totals[3]+" unavailable_tips="+totals[4]);
	}

	/** Reads unique sample identities and filename-genus grouping without consulting labels. */
	private static ArrayList<Sample> readSamples(final String path, final int expected){
		final ArrayList<Sample> samples=new ArrayList<Sample>();
		final HashSet<String> ids=new HashSet<String>();
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 fields=new LineParser1('\t');
		boolean error;
		try{
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				fields.set(line);
				if(fields.terms()!=12){throw new IllegalArgumentException("Expected the12-column frozen samples manifest.");}
				final String id=fields.parseString(0), group=fields.parseString(3);
				if(!id.matches("g[0-9]{4}") || !ids.add(id) || group.isEmpty()){
					throw new IllegalArgumentException("Invalid or duplicate sample identity: "+id);
				}
				samples.add(new Sample(id, group));
			}
		}finally{error=file.close();}
		if(error || samples.size()!=expected){throw new IllegalArgumentException("Incomplete samples manifest: "+samples.size());}
		return samples;
	}

	/** Fixed shuffled80/10/10 group allocation prevents rows from one genus crossing splits. */
	private static void assignSplits(final ArrayList<Sample> samples, final long seed){
		assert(!samples.isEmpty()) : "Split assignment requires at least one verified sample.";
		final HashSet<String> unique=new HashSet<String>();
		for(Sample sample : samples){unique.add(sample.group);}
		final ArrayList<String> groups=new ArrayList<String>(unique);
		Collections.sort(groups);
		Collections.shuffle(groups, new Random(seed));
		final HashMap<String, Integer> assignments=new HashMap<String, Integer>();
		for(int i=0; i<groups.size(); i++){
			assignments.put(groups.get(i), i<groups.size()*8/10 ? 0 : i<groups.size()*9/10 ? 1 : 2);
		}
		for(Sample sample : samples){sample.split=assignments.get(sample.group);}
	}

	/** Lockstep identity checks bind each original vector to its trace and unchanged reference label. */
	private static long[] convert(final Sample sample, final String root, final String out,
			final FusionTipEntropy entropy, final ByteStreamWriter training){
		assert(sample.split>=0 && sample.split<3) : "Each genome must have a label-independent split before conversion.";
		final String dir=root+"/data/"+sample.id;
		final ByteFile[] files=new ByteFile[3];
		ByteStreamWriter vectors=null;
		final long[] counts=new long[5];
		boolean error=false;
		try{
			files[0]=ByteFile.makeByteFile(dir+"/joins.vectors.tsv", true);
			files[1]=ByteFile.makeByteFile(dir+"/joins.trace.tsv", true);
			files[2]=ByteFile.makeByteFile(dir+"/labels.tsv", true);
			requireLine(files[0], "#schema="+FusionJoinFeatures.VERSION);
			final ByteBuilder header=new ByteBuilder("event\tphase_k\tquery_k");
			for(String name : FusionJoinFeatures.names()){header.tab().append(name);}
			requireLine(files[0], header.toString());
			requireLine(files[2], "#counts_saturated_at=2; target=NA is excluded from binary training");
			requireLine(files[2], "event\tphase_k\tquery_k\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits");
			vectors=open(out+"/"+sample.id+".vectors.tsv");
			vectors.print("#schema="+FusionTipEntropy.VERSION+"\n");
			vectors.print(header.tab().append(FusionTipEntropy.COLUMNS).nl());
			final LineParser1 v=new LineParser1('\t'), trace=new LineParser1('\t'), label=new LineParser1('\t');
			final ByteBuilder row=new ByteBuilder(), extra=new ByteBuilder();
			int phase=-1;
			for(byte[] line=files[0].nextLine(); line!=null; line=files[0].nextLine()){
				v.set(line);
				if(v.terms()!=110 || v.parseInt(0)!=counts[0]+1 || v.parseInt(1)!=v.parseInt(2)){
					throw new IllegalArgumentException("Invalid same-K vector identity in "+sample.id);
				}
				for(int i=3; i<110; i++){
					if(!Float.isFinite(v.parseFloat(i))){throw new IllegalArgumentException("Nonfinite source vector.");}
				}
				boolean found=false;
				for(byte[] t=files[1].nextLine(); t!=null; t=files[1].nextLine()){
					if(t.length==0 || t[0]=='#'){continue;}
					trace.set(t);
					if(trace.termEquals("FusionPass", 0)){phase=trace.parseInt(1);}
					else if(trace.termEquals("OverlapPair", 0)){found=true; break;}
					else{throw new IllegalArgumentException("Unknown trace record.");}
				}
				if(!found || trace.terms()!=11 || phase!=v.parseInt(1)){
					throw new IllegalArgumentException("Missing or mismatched trace for "+sample.id+" event "+v.parseInt(0));
				}
				final byte[] target=files[2].nextLine();
				if(target==null){throw new IllegalArgumentException("Missing label in "+sample.id);}
				label.set(target);
				if(label.terms()!=8 || label.parseInt(0)!=v.parseInt(0) ||
						label.parseInt(1)!=phase || label.parseInt(2)!=phase){
					throw new IllegalArgumentException("Label identity mismatch in "+sample.id);
				}
				extra.clear();
				counts[4]+=entropy.append(extra, trace.parseByteArray(9), trace.parseByteArray(10));
				vectors.print(row.clear().append(line).append(extra).nl());
				counts[0]++;
				final int value;
				if(label.termEquals("supported", 3) && label.termEquals("1", 4)){value=1; counts[1]++;}
				else if(label.termEquals("contradicted", 3) && label.termEquals("0", 4)){value=0; counts[2]++;}
				else if(label.termEquals("unresolved", 3) && label.termEquals("NA", 4)){counts[3]++; continue;}
				else{throw new IllegalArgumentException("Invalid truth label in "+sample.id);}
				int start=0;
				for(int tabs=0; tabs<3; start++){if(line[start]=='\t'){tabs++;}}
				training.print(row.clear().append(line, start, line.length-start).append(extra).tab().append(value).nl());
			}
			if(files[2].nextLine()!=null){throw new IllegalArgumentException("Extra labels in "+sample.id);}
			for(byte[] line=files[1].nextLine(); line!=null; line=files[1].nextLine()){
				trace.set(line);
				if(trace.termEquals("OverlapPair", 0)){throw new IllegalArgumentException("Extra trace pairs in "+sample.id);}
			}
		}finally{
			for(ByteFile file : files){if(file!=null){error=file.close() | error;}}
			if(vectors!=null){error=vectors.poisonAndWait() | error;}
		}
		if(error){throw new IllegalStateException("Incomplete entropy output for "+sample.id);}
		return counts;
	}

	/** Requires exact versioned headers; a similarly shaped future schema must not silently pass. */
	private static void requireLine(final ByteFile file, final String expected){
		final byte[] line=file.nextLine();
		if(line==null || !expected.equals(new String(line, StandardCharsets.UTF_8))){
			throw new IllegalArgumentException("Unexpected corpus header; expected "+expected);
		}
	}

	/** Opens only new regular output paths inside the fresh derivative directory. */
	private static ByteStreamWriter open(final String path){
		if(!Tools.testOutputFiles(false, false, false, path)){throw new IllegalArgumentException("Output exists: "+path);}
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/** Immutable identity and assigned split; no reference label participates in assignment. */
	private static final class Sample {
		/** Retains only the manifest fields needed for deterministic grouping and path resolution. */
		Sample(final String id_, final String group_){id=id_; group=group_;}
		final String id, group;
		int split=-1;
	}

	private static final String[] SPLITS={"train", "validate", "test"};
}
