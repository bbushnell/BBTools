package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Private historical-mask evidence preparation. Original annotations and truth stay immutable.
 * Each target is one genomic-forward window; no strand reversal or inferred coordinates.
 * @author Brian, Raiden */
public final class EukPanelHistoricalReview {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelHistoricalReview.class, false);
		String masks=null, inputs=null, sources=null, contigs=null, contigSources=null, models=null, out=null;
		int expected=724, assemblies=198, flank=500;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("masks")){masks=v;}else if(k.equals("inputs")){inputs=v;}else if(k.equals("sources")){sources=v;}
			else if(k.equals("contigs")){contigs=v;}else if(k.equals("contigsources")){contigSources=v;}else if(k.equals("models")){models=v;}
			else if(k.equals("outdir")){out=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else if(k.equals("assemblies")){assemblies=Integer.parseInt(v);}else if(k.equals("flank")){flank=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(masks!=null && inputs!=null && sources!=null && contigs!=null && contigSources!=null && models!=null && out!=null && expected>0 && assemblies>0 && flank>=0, "Explicit bound sources, counts, flank and fresh outdir required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelHistoricalReview x=new EukPanelHistoricalReview();
		x.panel.hashes=EukPanelGffValidation.hashes(sources);x.panel.contigHashes=EukPanelGffValidation.hashes(contigSources);
		x.panel.loadInputs(inputs, assemblies);x.panel.loadContigs(contigs);x.loadInputs(inputs);x.loadMasks(masks, expected, flank);
		require(new File(out).mkdir(), "Fresh outdir required: "+out);
		x.models(models, out);x.extract(out);x.write(out);Shared.closeStream(pp.outstream);
	}
	void loadInputs(String path){
		assert(fastas.isEmpty()):"One assembly-qualified FASTA binding per panel member";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{in.nextLine();for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==15, "Input schema changed");require(fastas.put(p.parseString(1), p.parseString(6))==null, "Duplicate assembly");}}
		finally{require(!in.close(), "Input read failed");}
	}
	void loadMasks(String path, int expected, int flank){
		assert(rows.isEmpty()):"Retain exactly the declared eligible-mask population";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in, EukPanelHistoricalMask.OUTPUT_HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==14 && (p.termEquals("0", 7) || p.termEquals("1", 7)), "Mask schema/eligibility changed");
			if(p.termEquals("0", 7)){continue;}
			final String assembly=p.parseString(0), seqid=p.parseString(1), strand=p.parseString(2), source=p.parseString(8);
			require(seen.add(assembly+'\t'+source+'\t'+p.parseLong(10)), "Duplicate assembly-qualified source identity");
			require(p.parseString(9).equals(panel.hashes.get(source)), "Original annotation source hash must match manifest");
			require(p.termEquals("HIST_ONLY", 6) || p.termEquals("HIST_LINKED_UNSCORED", 6), "Admitted annotations are outside this review");
			final EukPanelGffValidation.Assembly a=panel.byAssembly.get(assembly);require(a!=null, "Unknown assembly");
			final EukPanelGffValidation.Contig c=a.literal.get(seqid);require(c!=null, "Eligible mask must resolve to its exact native contig");
			final long start=p.parseLong(3), end=p.parseLong(4);require(start>=1 && end>=start && end<=c.length && (strand.equals("+") || strand.equals("-")), "Eligible mask interval/strand is invalid; do not repair silently");
			final Row r=new Row("H"+(rows.size()+1), assembly, seqid, strand, p.parseString(5), line.clone(), start, end, c.length, c.compartment, fastas.get(assembly), flank);
			rows.add(r);families.add(r.family);
			ArrayList<Row> list=byAssembly.get(assembly);if(list==null){list=new ArrayList<Row>();byAssembly.put(assembly, list);}list.add(r);
		}}finally{require(!in.close(), "Mask read failed");}
		require(rows.size()==expected, "Eligible mask count changed: "+rows.size()+" expected="+expected);
	}
	/** Split complete CM+filter-HMM records; preserve every byte except newline normalization. */
	void models(String path, String out){
		assert(!families.isEmpty()):"Only families present in the declared masks receive scan bundles";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectMap<String,String> accessionFamily=new ObjectMap<String,String>(String.class, String.class);
		final ArrayList<String> bundles=new ArrayList<String>();final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final ObjectMap<String,ByteStreamWriter> writers=new ObjectMap<String,ByteStreamWriter>(String.class, ByteStreamWriter.class);
		final ObjectSet<String> written=new ObjectSet<String>(String.class);final ArrayList<String> expected=new ArrayList<String>();
		try{RrnaResourceIO.header(in, MODEL_HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==8, "Model schema changed");
			final String acc=p.parseString(1), family=p.parseString(2), source=p.parseString(6);
			require(accessionFamily.put(acc, family)==null, "Repeated accession requires explicit byte arbitration");
			if(families.contains(family)){expected.add(acc);if(paths.add(source)){bundles.add(source);}}
		}}finally{require(!in.close(), "Model metadata read failed");}
		try{
			for(String family:families.keys()){if(family!=null){
				boolean found=false;for(String acc:expected){found|=family.equals(accessionFamily.get(acc));}
				require(found, "No independent model for mask family: "+family);
				writers.put(family, RrnaResourceIO.writer(out+'/'+family+".cm"));
			}}
			for(String pathCM:bundles){
			final ByteFile cm=RrnaResourceIO.open(pathCM);final ByteBuilder record=new ByteBuilder();String acc=null;boolean header=false;
			try{for(byte[] line;(line=cm.nextLine())!=null;){
				if(EukOrganellePanelInventory.starts(line, "INFERNAL1/")){
					if(record.length()>0){emitModel(record, acc, accessionFamily, writers, written);}record.clear();acc=null;header=true;
				}
				if(header && EukOrganellePanelInventory.starts(line, "ACC ")){require(acc==null, "Duplicate CM accession");acc=EukOrganellePanelInventory.word(line, 1);}
				if(header && EukOrganellePanelInventory.starts(line, "CM")){header=false;}
				record.append(line).nl();
			}if(record.length()>0){emitModel(record, acc, accessionFamily, writers, written);}}
			finally{require(!cm.close(), "Model bundle read failed");}
		}for(String acc:expected){require(written.contains(acc), "Expected model absent from complete-record split: "+acc);}}
		finally{for(ByteStreamWriter w:writers.values()){if(w!=null){RrnaResourceIO.finish(w);}}}
		for(String family:families.keys()){if(family!=null){require(writers.contains(family), "No independent model for mask family: "+family);}}
	}
	static void emitModel(ByteBuilder record, String acc, ObjectMap<String,String> families, ObjectMap<String,ByteStreamWriter> writers, ObjectSet<String> seen){
		require(acc!=null && families.contains(acc), "Unrecognized or incomplete CM record");
		final ByteStreamWriter w=writers.get(families.get(acc));if(w!=null){require(seen.add(acc), "Duplicate selected CM");w.print(record);}
	}
	void extract(String out){
		assert(!rows.isEmpty()):"No implicit whole-genome scans; extract only declared local windows";
		for(EukPanelGffValidation.Assembly a:panel.assemblies){
			final ArrayList<Row> list=byAssembly.get(a.id);if(list==null){continue;}
			final ObjectMap<String,ArrayList<Row>> contigs=new ObjectMap<String,ArrayList<Row>>(String.class, (Class)ArrayList.class);
			for(Row r:list){ArrayList<Row> group=contigs.get(r.seqid);if(group==null){group=new ArrayList<Row>();contigs.put(r.seqid, group);}group.add(r);}
			RrnaResourceIO.stream(fastas.get(a.id), read->{
				final ArrayList<Row> group=contigs.get(RrnaResourceIO.first(read.id));if(group==null){return;}
				for(Row r:group){require(!r.written && read.length()==r.contigLength, "FASTA identity/length differs from bound contig evidence");
					final ByteStreamWriter w=RrnaResourceIO.writer(out+'/'+r.id+".fa");final ByteBuilder b=new ByteBuilder();
					try{w.print(b.append('>').append(r.id).nl());for(long i=r.left-1;i<r.right;i+=80){b.clear();final int len=(int)Math.min(80, r.right-i);w.print(b.append(read.bases, (int)i, len).nl());}}
					finally{RrnaResourceIO.finish(w);}r.written=true;
				}
			});
		}
		for(Row r:rows){require(r.written, "Requested source window not found: "+r.id);}
	}
	void write(String out){
		assert(rows.size()>0):"Every completed window needs an assembly/source-qualified mapping";
		final ByteStreamWriter table=RrnaResourceIO.writer(out+"/windows.tsv"), bindings=RrnaResourceIO.writer(out+"/source_files.sha80"), tasks=RrnaResourceIO.writer(out+"/tasks.tsv");
		final ByteBuilder b=new ByteBuilder();final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{table.println(WINDOW_HEADER);tasks.println("id\tfamily\tfasta\tmodels\twindow_length");
			for(Row r:rows){table.print(b.clear().append(r.id).tab().append(r.original).tab().append(r.contigLength).tab().append(r.compartment).tab().append(r.left).tab().append(r.right).tab().append(r.fasta).nl());
				tasks.print(b.clear().append(r.id).tab().append(r.family).tab().append(out+'/'+r.id+".fa").tab().append(out+'/'+r.family+".cm").tab().append(r.right-r.left+1).nl());
				if(seen.add(r.fasta)){bindings.print(b.clear().append(panel.hashes.get(r.fasta)).tab().append(r.fasta).nl());}
			}
		}finally{RrnaResourceIO.finish(table, bindings, tasks);}
		System.err.println("HISTORICAL_WINDOWS_COMPLETE masks="+rows.size()+" source_fastas="+seen.size()+" frozen_truth_changed=false");
	}
	public static long left(long start, int flank){require(start>=1 && flank>=0, "One-based source start and nonnegative flank required");return Math.max(1, start-flank);}
	public static long right(long end, long length, int flank){require(end>=1 && end<=length && flank>=0, "Valid source end required");return end+Math.min(length-end, flank);}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Row{
		Row(String i,String a,String s,String st,String f,byte[] raw,long from,long to,long len,String comp,String fa,int flank){id=i;assembly=a;seqid=s;strand=st;family=f;original=raw;start=from;end=to;contigLength=len;compartment=comp;fasta=fa;left=left(from,flank);right=right(to,len,flank);}
		final String id,assembly,seqid,strand,family,compartment,fasta;final byte[] original;final long start,end,contigLength,left,right;boolean written;
	}
	final EukPanelGffValidation panel=new EukPanelGffValidation();
	final ObjectMap<String,String> fastas=new ObjectMap<String,String>(String.class,String.class);
	final ObjectMap<String,ArrayList<Row>> byAssembly=new ObjectMap<String,ArrayList<Row>>(String.class,(Class)ArrayList.class);
	final ArrayList<Row> rows=new ArrayList<Row>();final ObjectSet<String> families=new ObjectSet<String>(String.class);
	public static final String MODEL_HEADER="bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80";
	public static final String WINDOW_HEADER="mask_id\t"+EukPanelHistoricalMask.OUTPUT_HEADER+"\tcontig_length\tcompartment\twindow_start\twindow_end\tsource_fasta";
}
