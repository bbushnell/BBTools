package prok;

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

/** Validates historical intervals without altering, dropping or admitting source features.
 * @author Raiden */
public final class EukPanelGffValidation {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelGffValidation.class, false);
		String inputs=null, sources=null, contigs=null, contigSources=null, out=null;int expected=198;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("inputs")){inputs=v;}else if(k.equals("sources")){sources=v;}else if(k.equals("contigs")){contigs=v;}
			else if(k.equals("contigsources")){contigSources=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}
		}
		require(inputs!=null && sources!=null && contigs!=null && contigSources!=null && out!=null && expected>0,
			"Bind corrected annotation inputs and observed-contig source hashes before validation");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);ReadWrite.ZIPLEVEL=1;
		final EukPanelGffValidation tool=new EukPanelGffValidation();
		tool.hashes=hashes(sources);tool.contigHashes=hashes(contigSources);
		tool.loadInputs(inputs, expected);tool.loadContigs(contigs);tool.write(out);Shared.closeStream(pp.outstream);
	}
	public static ObjectMap<String,String> hashes(String path){
		assert(path!=null):"Every source file is tied to a verified provenance manifest";
		final ObjectMap<String,String> map=new ObjectMap<String,String>(String.class, String.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] line; (line=in.nextLine())!=null; ){
			p.set(line);require(p.terms()==2, "Source binding must be sha80/path");
			final String hash=p.parseString(0), file=p.parseString(1);
			require(SeedOffsetTable.validSha80(hash) && map.put(file, hash)==null, "Source hashes must be valid and paths unique");
		}}finally{require(!in.close(), "Source hash manifest read failed");}
		return map;
	}
	public void loadInputs(String path, int expected){
		assert(assemblies.isEmpty()):"Load exactly one approved input set";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			p.set(in.nextLine());require(p.terms()==15 && p.termEquals("assembly", 1) && p.termEquals("gff", 8), "Expected corrected input binding schema");
			for(byte[] line; (line=in.nextLine())!=null; ){
				p.set(line);require(p.terms()==15, "Input binding schema changed");
				final String fasta=p.parseString(6), digest=hashes.get(fasta), gff=p.parseString(8);
				require(digest!=null && digest.equals(contigHashes.get(fasta)), "Annotation binding and contig observation must use identical FASTA bytes: "+fasta);
				final String selected=p.parseString(7);require(selected.equals("NA") || selected.equals(digest), "Frozen selected FASTA hash differs");
				final Assembly a=new Assembly(p.parseString(0), p.parseString(1), gff, p.parseLong(4), p.parseLong(5));
				a.gffHash=gff.equals("NA")?"NA":hashes.get(gff);
				require(a.gffHash!=null && a.bases>0 && a.records>0 && byAssembly.put(a.id, a)==null, "Unique assemblies need positive dimensions and bound annotation files");
				assemblies.add(a);
			}
		}finally{require(!in.close(), "Corrected input manifest read failed");}
		require(assemblies.size()==expected, "Approved assembly count changed: "+assemblies.size()+" expected="+expected);
	}
	void loadContigs(String path){
		assert(!assemblies.isEmpty()):"Observed contigs cannot invent assemblies outside the selected panel";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final int[] tid=new int[1];
		try{
			p.set(in.nextLine());require(p.terms()==15 && p.termEquals("seqid", 2) && p.termEquals("native_taxid", 12), "Expected evidence-backed contig schema");
			for(byte[] line; (line=in.nextLine())!=null; ){
				p.set(line);require(p.terms()==15, "Contig observation schema changed");
				final Assembly a=byAssembly.get(p.parseString(1));require(a!=null && p.termEquals(a.group, 0), "Contig assembly/group absent from bound inputs");
				final String id=p.parseString(2), bare=EukPanelCompartments.accession(id, tid), evidence=p.parseString(5), reported=p.parseString(8);
				final Contig c=new Contig(id, p.parseLong(3), p.parseString(4), tid[0]);
				require(c.length>0 && tid[0]==p.parseInt(12) && a.literal.put(id, c)==null, "Contig IDs, native tags and positive observed lengths must be unambiguous");
				a.observedRecords++;a.observedBases+=c.length;alias(a, bare, c, "VALIDATED_NATIVE_TAG_ACCESSION");
				if(evidence.equals("ASSEMBLY_REPORT_GENBANK_EQUAL") && !reported.equals("NA")){
					alias(a, reported, c, "EXACT_REPORT_EQUAL_ALIAS");
				}
			}
		}finally{require(!in.close(), "Observed contig table read failed");}
		for(Assembly a:assemblies){require(a.observedRecords==a.records && a.observedBases==a.bases,
			"Contig table must conserve the complete bound assembly: "+a.id);}
	}
	static void alias(Assembly a, String name, Contig c, String kind){
		assert(name!=null && c!=null):"Only explicit observed-contig evidence creates an alias";
		final Alias old=a.aliases.get(name);
		if(old==null){a.aliases.put(name, new Alias(c, kind));}
		else if(old.contig!=c){a.ambiguous.add(name);}
	}
	static Resolution resolve(Assembly a, String seqid){
		assert(a!=null && seqid!=null):"Resolve within one assembly; accessions alone are not globally unique";
		final int[] tid=new int[1];final String bare;
		try{bare=EukPanelCompartments.accession(seqid, tid);}catch(IllegalArgumentException malformed){return new Resolution(null, "MALFORMED_NATIVE_TAG", 0);}
		final Contig literal=a.literal.get(seqid);
		if(literal!=null){return new Resolution(literal, "EXACT_LITERAL", tid[0]);}
		if(a.ambiguous.contains(bare)){return new Resolution(null, "AMBIGUOUS_SEQID", tid[0]);}
		final Alias alias=a.aliases.get(bare);
		return alias==null?new Resolution(null, "UNRESOLVED_SEQID", tid[0]):new Resolution(alias.contig, alias.kind, tid[0]);
	}
	/** Invalid coordinates remain source evidence; they never become invented zeroes. */
	static long coordinate(LineParser1 p, int term){
		p.setBounds(term);final byte[] line=p.line();final int start=p.a(), end=p.b();
		if(start==end){return Long.MIN_VALUE;}
		long n=0;int i=start;final boolean negative=line[i]=='-';if(negative){i++;}
		if(i==end){return Long.MIN_VALUE;}
		for(; i<end; i++){
			final int digit=line[i]-'0';
			if(digit<0 || digit>9 || n>(Long.MAX_VALUE-digit)/10){return Long.MIN_VALUE;}
			n=n*10+digit;
		}
		return negative?-n:n;
	}
	void write(String prefix){
		assert(!assemblies.isEmpty()):"Historical validation retains every selected assembly, including unavailable annotations";
		RrnaResourceIO.fresh(prefix, ".features.tsv.gz", ".summary.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".features.tsv.gz");ByteStreamWriter summary=null;
		final ByteBuilder b=new ByteBuilder();long total=0;
		try{
			summary=RrnaResourceIO.writer(prefix+".summary.tsv");
			out.println("selection_group\tassembly\tsource_gff\tsource_sha80\tsource_line\tvalidation_status\tmapping_status\trange_status\tnative_seqid\tcontig_length\tcompartment\tgff_taxid\tnative_taxid\ttaxonomy_status\toriginal_seqid\toriginal_source\toriginal_type\toriginal_start\toriginal_end\toriginal_score\toriginal_strand\toriginal_phase\toriginal_attributes");
			summary.println("selection_group\tassembly\tsource_gff\tsource_sha80\tavailability\tsource_lines\tcomments\tblank_lines\tinline_fasta_lines\tfeatures\tvalid_intervals\tunresolved_seqids\tout_of_range\tinvalid_coordinates");
			for(Assembly a:assemblies){
				final long[] counts=new long[4];long lines=0, comments=0, blanks=0, embedded=0, features=0;
				if(!a.gff.equals("NA")){
					final ByteFile in=RrnaResourceIO.open(a.gff);final LineParser1 p=new LineParser1('\t');
					boolean inline=false;String last=null;Resolution r=null;
					try{for(byte[] line; (line=in.nextLine())!=null; ){
						lines++;
						if(inline){embedded++;continue;}
						if(line.length==0){blanks++;continue;}
						if(line[0]=='#'){
							comments++;if(EukOrganellePanelInventory.starts(line, "##FASTA") && (line.length==7 || line[7]<=32)){inline=true;}continue;
						}
						p.set(line);require(p.terms()==9, "Historical GFF feature requires nine fields; source retained at "+a.gff+":"+lines);
						if(last==null || !p.termEquals(last, 0)){last=p.parseString(0);r=resolve(a, last);}
						final long start=coordinate(p, 3), stop=coordinate(p, 4);
						final boolean invalid=start<1 || stop<start;
						final String range=invalid?"INVALID_COORDINATES":r.contig==null?"UNKNOWN_CONTIG_LENGTH":stop>r.contig.length?"OUT_OF_RANGE":"VALID_1_BASED_INCLUSIVE";
						final int state=invalid?3:r.contig==null?1:stop>r.contig.length?2:0;
						counts[state]++;features++;
						b.clear().append(a.group).tab().append(a.id).tab().append(a.gff).tab().append(a.gffHash).tab().append(lines).tab().append(STATES[state]).tab().append(r.kind).tab().append(range);
						if(r.contig==null){b.append("\tNA\tNA\tUNKNOWN").tab().append(r.taxid).append("\t0\tUNRESOLVED");}
						else{
							final Contig c=r.contig;final String taxonomy=r.taxid==0?"GFF_NOT_TAGGED":c.taxid==0?"FASTA_NOT_TAGGED":r.taxid==c.taxid?"MATCH":"MISMATCH";
							b.tab().append(c.id).tab().append(c.length).tab().append(c.compartment).tab().append(r.taxid).tab().append(c.taxid).tab().append(taxonomy);
						}
						out.print(b.tab().append(line).nl()); // All nine original fields and attributes remain byte-for-byte intact.
					}}finally{require(!in.close(), "Historical GFF read failed: "+a.gff);}
				}
				require(lines==comments+blanks+embedded+features && features==counts[0]+counts[1]+counts[2]+counts[3], "Every source line and feature must be conserved: "+a.id);
				b.clear().append(a.group).tab().append(a.id).tab().append(a.gff).tab().append(a.gffHash).tab().append(a.gff.equals("NA")?"MISSING_GFF":"BOUND_SOURCE").tab().append(lines).tab().append(comments).tab().append(blanks).tab().append(embedded).tab().append(features);
				for(long count:counts){b.tab().append(count);}summary.print(b.nl());total+=features;
			}
		}finally{if(summary==null){RrnaResourceIO.finish(out);}else{RrnaResourceIO.finish(out, summary);}}
		System.err.println("HISTORICAL_GFF_VALIDATION_COMPLETE assemblies="+assemblies.size()+" source_features="+total+" no_truth_admission=true no_coordinate_repair=true");
	}
	public static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final String[] STATES={"VALID_INTERVAL", "UNRESOLVED_SEQID", "OUT_OF_RANGE", "INVALID_COORDINATES"};
	static final class Contig{Contig(String i, long l, String c, int t){id=i;length=l;compartment=c;taxid=t;}final String id,compartment;final long length;final int taxid;}
	static final class Alias{Alias(Contig c, String k){contig=c;kind=k;}final Contig contig;final String kind;}
	static final class Resolution{Resolution(Contig c, String k, int t){contig=c;kind=k;taxid=t;}final Contig contig;final String kind;final int taxid;}
	static final class Assembly{
		Assembly(String g, String i, String f, long b, long r){group=g;id=i;gff=f;bases=b;records=r;}
		final String group,id,gff;final long bases,records;String gffHash;long observedBases,observedRecords;
		final ObjectMap<String,Contig> literal=new ObjectMap<String,Contig>(String.class, Contig.class);
		final ObjectMap<String,Alias> aliases=new ObjectMap<String,Alias>(String.class, Alias.class);
		final ObjectSet<String> ambiguous=new ObjectSet<String>(String.class);
	}
	public ObjectMap<String,String> hashes,contigHashes;
	final ObjectMap<String,Assembly> byAssembly=new ObjectMap<String,Assembly>(String.class, Assembly.class);
	public final ArrayList<Assembly> assemblies=new ArrayList<Assembly>();
}
