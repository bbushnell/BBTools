package prot;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import clade.CladeSearcher;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import tax.TaxNode;
import tax.TaxTree;
import template.ThreadWaiter;

/**
 * Compares the assembly client's taxonomy path against the frozen R1 replay.
 * Reads the original FASTAs without gene calling or MAG-QC quality-network inference.
 * Normal QuickClade ranking and confidence networks remain part of taxonomy.
 * @author Brian Bushnell, Yoimiya
 */
public final class MagQCNormalTaxonomyParity {

	/** Requires explicit input pins and either local classification or a private loopback server. */
	public static void main(String[] args){
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg:new PreParser(args, null, false).args){
			final int eq=arg.indexOf('=');
			if(eq<1 || options.put(arg.substring(0, eq).toLowerCase(Locale.ROOT), arg.substring(eq+1))!=null){
				throw new IllegalArgumentException("Malformed or duplicate option: "+arg);
			}
		}
		final String truth=take(options, "truth"), truthPin=take(options, "truthsha80");
		final String expected=take(options, "expected"), expectedPin=take(options, "expectedsha80");
		final String treePath=take(options, "tree"), treePin=take(options, "treesha80");
		final String directory=take(options, "fastadir"), output=take(options, "out");
		final int threads=Integer.parseInt(take(options, "t"));
		if(threads<1 || threads>32){throw new IllegalArgumentException("t must be 1..32");}
		Shared.setThreads(threads);
		final HashMap<String,String> taxonomy=new HashMap<String,String>();
		taxonomy.put("taxmode", take(options, "taxmode"));
		taxonomy.put("normalsearch", take(options, "normalsearch"));
		final boolean normal=MagQCAssemblyInput.normalSearch(taxonomy);
		if(MagQCAssemblyInput.localTaxonomy(taxonomy)){
			taxonomy.put("taxsketch", take(options, "taxsketch"));
			pin(taxonomy.get("taxsketch"), take(options, "taxsketchsha80"));
			final String ref=CladeSearcher.defaultRef();
			if(ref==null){throw new IllegalArgumentException("No local QuickClade reference");}
			pin(ref, take(options, "refsha80"));
		}else if("server".equals(taxonomy.get("taxmode"))){
			final String address=take(options, "taxaddress");
			if(!address.matches("http://127\\.0\\.0\\.1:[0-9]+(?:/[^\\s]*)?")){
				throw new IllegalArgumentException("Parity tests require a private loopback server: "+address);
			}
			taxonomy.put("taxaddress", address);
		}else{throw new IllegalArgumentException("taxmode must be local or server");}
		if(!options.isEmpty()){throw new IllegalArgumentException("Unexpected options: "+options.keySet());}
		pin(truth, truthPin); pin(expected, expectedPin); pin(treePath, treePin);
		final ArrayList<Bin> bins=loadInputs(truth, expected, directory, normal);
		try(MagQCAssemblyInput.SketchSession session=MagQCAssemblyInput.openSketchSession(normal)){
			final ArrayList<Worker> workers=new ArrayList<Worker>();
			for(int i=0; i<threads; i++){workers.add(new Worker(bins, session, taxonomy, i, threads));}
			if(!ThreadWaiter.startAndWait(workers)){throw new IllegalStateException("Taxonomy workers did not join");}
			for(Worker worker:workers){
				if(worker.failure!=null){throw new RuntimeException("Product taxonomy parity worker failed", worker.failure);}
			}
		}
		// Loading a TaxTree publishes shared state read by CladeIndex.sketchLCA.
		// Grading must not change the candidate features used by the product call.
		final TaxTree tree=TaxTree.loadTaxTree(treePath, System.err, false, false);
		final int[] counts=new int[3];
		int mismatches=0;
		final ByteBuilder rows=new ByteBuilder("bin\texpected_reference_tid\tactual_reference_tid\tprimary_phylum_tid\tactual_phylum_tid\tgroup\tproduct_domain\tproduct_phylum\n");
		for(Bin bin:bins){
			if(bin.actual==null){throw new IllegalStateException("Missing product result: "+bin.id);}
			final int phylum=phylum(tree, bin.actual.taxId);
			final int group=bin.primary<1 || phylum<1 ? 2 : phylum==bin.primary ? 0 : 1;
			counts[group]++;
			if(bin.actual.taxId!=bin.expected || phylum!=bin.expectedPhylum){mismatches++;}
			rows.append(bin.id).tab().append(bin.expected).tab().append(bin.actual.taxId)
				.tab().append(bin.primary).tab().append(phylum).tab().append(GROUPS[group])
				.tab().append(bin.actual.domain).tab().append(bin.actual.phylum).nl();
		}
		write(output, rows);
		final int[] wanted=normal ? new int[]{932, 55, 13} : new int[]{742, 233, 25};
		if(mismatches!=0 || !java.util.Arrays.equals(counts, wanted)){
			throw new IllegalStateException("Taxonomy parity differs: mismatches="+mismatches+
				" counts="+java.util.Arrays.toString(counts)+"; inspect "+output);
		}
		System.out.println("NORMAL_TAXONOMY_PARITY_PASS bins=1000 transport="+taxonomy.get("taxmode")+
			" normalsearch="+normal+" correct="+counts[0]+" incorrect="+counts[1]+" unresolved="+counts[2]);
	}

	/** Joins immutable FASTA metadata to the independently measured replay references. */
	private static ArrayList<Bin> loadInputs(String truth, String expected, String directory, boolean normal){
		final ArrayList<Bin> bins=new ArrayList<Bin>();
		final ObjectIntMap<String> indexes=new ObjectIntMap<String>(String.class);
		try(Table table=new Table(truth)){
			while(table.next()){
				if(!table.text("regime").equals("R1")){continue;}
				final String id=table.text("bin_id");
				if(!id.matches("bp_R1_[0-9]{4}") || indexes.contains(id)){throw new IllegalArgumentException("Unexpected/duplicate R1 bin: "+id);}
				indexes.put(id, bins.size());
				bins.add(new Bin(id, new File(directory, id+".fa").toString(), table.text("fasta_sha80"), table.number("bin_bp")));
			}
		}
		if(bins.size()!=1000){throw new IllegalArgumentException("Expected exactly 1000 R1 FASTAs");}
		try(Table table=new Table(expected)){
			while(table.next()){
				if(!table.text("mode").equals(normal ? "normal_local" : "panel_faithful")){continue;}
				final String id=table.text("bin"); final int index=indexes.get(id);
				if(index<0 || bins.get(index).seen){throw new IllegalArgumentException("Unexpected/duplicate expected bin: "+id);}
				final Bin bin=bins.get(index); bin.seen=true;
				bin.expected=table.number("replay_reference_tid");
				bin.primary=table.number("primary_phylum_tid");
				bin.expectedPhylum=table.number("replay_phylum_tid");
			}
		}
		for(Bin bin:bins){if(!bin.seen){throw new IllegalArgumentException("Missing expected bin: "+bin.id);}}
		return bins;
	}

	/** Each worker owns assembly reads; classification uses the product's immutable session policy. */
	private static final class Worker extends Thread{
		Worker(ArrayList<Bin> bins_, MagQCAssemblyInput.SketchSession session_, HashMap<String,String> options_, int first_, int stride_){
			bins=bins_; session=session_; options=options_; first=first_; stride=stride_;
			assert(stride>0 && first>=0 && first<stride) : "Worker strides must partition the fixed input list";
		}
		@Override public void run(){
			try{
				for(int i=first; i<bins.size(); i+=stride){
					final Bin bin=bins.get(i); pin(bin.path, bin.pin);
					final ArrayList<Read> reads=MagQCAssemblyInput.readContigs(bin.path);
					long bases=0; for(Read read:reads){bases+=read.length();}
					if(bases!=bin.bases){throw new IllegalArgumentException("FASTA base count changed: "+bin.id);}
					bin.actual=session.classify(reads, options);
				}
			}catch(Throwable t){failure=t;}
		}
		private final ArrayList<Bin> bins;
		private final MagQCAssemblyInput.SketchSession session;
		private final HashMap<String,String> options;
		private final int first, stride;
		private Throwable failure;
	}

	/** Resolves exactly the same NCBI phylum level as the frozen replay assessment. */
	private static int phylum(TaxTree tree, long tid){
		if(tid<1 || tid>Integer.MAX_VALUE){return -1;}
		final TaxNode start=tree.getNode((int)tid, true);
		final TaxNode node=start==null ? null : tree.getNodeAtLevelExtended(start.id, TaxTree.PHYLUM_E);
		return node!=null && node.levelExtended==TaxTree.PHYLUM_E ? node.id : -1;
	}

	/** A small named-column reader; row parsing stays in reusable byte buffers. */
	private static final class Table implements AutoCloseable{
		Table(String path){
			file=ByteFile.makeByteFile(path, true);
			final byte[] header=file.nextLine();
			if(header==null){throw new IllegalArgumentException("Missing table header: "+path);}
			parser.set(header);
			for(int i=0; i<parser.terms(); i++){
				final String name=parser.parseString(i);
				if(columns.contains(name)){throw new IllegalArgumentException("Duplicate column: "+name);}
				columns.put(name, i);
			}
		}
		boolean next(){
			byte[] line;
			while((line=file.nextLine())!=null){
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				if(parser.terms()!=columns.size()){throw new IllegalArgumentException("Table width changed: "+file.name());}
				return true;
			}
			return false;
		}
		int column(String name){
			final int index=columns.get(name);
			if(index<0){throw new IllegalArgumentException("Missing column: "+name);}
			return index;
		}
		String text(String name){return parser.parseString(column(name));}
		long number(String name){return parser.parseLong(column(name));}
		@Override public void close(){if(file.close()){throw new IllegalStateException("Reader failed: "+file.name());}}
		private final ByteFile file;
		private final LineParser1 parser=new LineParser1('\t');
		private final ObjectIntMap<String> columns=new ObjectIntMap<String>(String.class);
	}

	/** Holds only durable input metadata and one immutable product result per bin. */
	private static final class Bin{
		Bin(String id_, String path_, String pin_, long bases_){id=id_; path=path_; pin=pin_; bases=bases_;}
		final String id, path, pin;
		final long bases;
		long primary, expected, expectedPhylum;
		boolean seen;
		MagQCAssemblyInput.Taxonomy actual;
	}

	/** Checks raw file identity on the cluster before accepting it as a frozen input. */
	private static void pin(String path, String expected){
		if(!expected.matches("[0-9a-f]{20}") || !DigestSuffix.file(path).equals(expected)){
			throw new IllegalArgumentException("Input pin mismatch: "+path);
		}
	}

	/** Removes a required argument so misspelled and unused options are rejected. */
	private static String take(HashMap<String,String> options, String key){
		final String value=options.remove(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing "+key);}
		return value;
	}

	/** Retains all mismatching rows for diagnosis and propagates writer failure. */
	private static void write(String path, ByteBuilder rows){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, true);
		writer.start(); writer.print(rows);
		if(writer.poisonAndWait()){throw new IllegalStateException("Output failed: "+path);}
	}

	private static final String[] GROUPS={"correct", "incorrect", "unresolved"};
}
