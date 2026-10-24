package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Assigns contig compartments from explicit report/GBFF evidence, never panel groups.
 * @author Raiden */
public final class EukPanelCompartments {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelCompartments.class, false);
		String inputs=null, sources=null, reports=null, reportHash=null, out=null;int expected=198;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("inputs")){inputs=value;}else if(key.equals("sources")){sources=value;}
			else if(key.equals("reports")){reports=value;}else if(key.equals("reportsha80")){reportHash=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("expected")){expected=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		require(inputs!=null && sources!=null && reports!=null && out!=null && SeedOffsetTable.validSha80(reportHash) && expected>0,
			"Explicit bound input table, source hashes, report stream/hash, fresh output and positive expected count are required");
		Shared.setThreads(1);
		final EukPanelCompartments mapper=new EukPanelCompartments();
		mapper.loadHashes(sources);mapper.loadInputs(inputs, expected);mapper.loadReports(reports, reportHash);
		for(String gbff:mapper.gbffPaths){mapper.loadGbff(gbff);}
		mapper.write(out);Shared.closeStream(pp.outstream);
	}

	void loadHashes(String path){
		assert(hashes.size()==0):"One already verified source manifest defines this run's files";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==2 && SeedOffsetTable.validSha80(p.parseString(0)), "Expected sha80/path source binding");
			final String digest=p.parseString(0), file=p.parseString(1);
			require(hashes.put(file, digest)==null, "Duplicate path in deduplicated source manifest: "+file);
		}}finally{require(!in.close(), "Source manifest read failed");}
	}
	void loadInputs(String path, int expected){
		assert(assemblies.isEmpty()):"The approved panel is loaded once, with no implicit additions";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> files=new ObjectSet<String>(String.class);
		try{
			p.set(in.nextLine());require(p.terms()==15 && p.termEquals("assembly", 1) && p.termEquals("fasta", 6)
				&& p.termEquals("gbff_record_accession", 11), "Expected reviewed input-binder schema");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==15, "Input row schema changed");
				final Assembly a=new Assembly(p.parseString(0), p.parseString(1), p.parseString(6), p.parseLong(4), p.parseLong(5), p.parseInt(3));
				require(a.group.equals("nuclear") || a.group.equals("plastid") || a.group.equals("mito"), "Unknown selection group");
				require(a.expectedBases>0 && a.expectedRecords>0 && hashes.contains(a.fasta), "FASTA requires positive dimensions and verified hash: "+a.id);
				final String selected=p.parseString(7);
				require(selected.equals("NA") || selected.equals(hashes.get(a.fasta)), "Selected and freshly bound FASTA hashes disagree: "+a.id);
				require(byAssembly.put(a.id, a)==null, "Ambiguous assembly identifier: "+a.id);assemblies.add(a);
				a.gbff=p.parseString(10);a.record=p.parseString(11);
				if(!a.gbff.equals("NA")){
					require(hashes.contains(a.gbff) && !a.record.equals("NA"), "Original GBFF/accession must have a verified source binding");
					require(byRecord.put(a.record, a)==null, "Duplicate original GBFF accession: "+a.record);
					if(files.add(a.gbff)){gbffPaths.add(a.gbff);}
				}
			}
		}finally{require(!in.close(), "Input binding read failed");}
		require(assemblies.size()==expected, "Conserve the approved input count; expected="+expected+" actual="+assemblies.size());
	}

	void loadReports(String path, String digest){
		assert(assemblies.size()>0):"Only reports for the already declared panel may contribute evidence";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');Assembly a=null;
		try{for(byte[] line;(line=in.nextLine())!=null;){
			if(EukOrganellePanelInventory.starts(line, "# Assembly name:")){a=null;}
			else if(EukOrganellePanelInventory.starts(line, "# RefSeq assembly accession:")){
				final String id=EukOrganellePanelInventory.value(line);a=byAssembly.get(id);
				if(a!=null){require(!a.reportSeen, "Duplicate selected assembly report: "+id);a.reportSeen=true;}
			}else if(a!=null && line.length>0 && line[0]!='#'){
				p.set(line);require(p.terms()==10, "NCBI report sequence row requires ten fields: "+a.id);
				final String ref=p.parseString(6), gb=p.parseString(4), relation=p.parseString(5), location=p.parseString(3), unit=p.parseString(7);
				final Evidence e=new Evidence(p.parseLong(8), reportCompartment(location, unit), "ASSEMBLY_REPORT", path, digest, ref, location, unit, "NA");
				require(e.length>0, "Report sequence length must be positive: "+a.id);
				if(!ref.equals("na")){putAlias(a, ref, e, "REFSEQ");}
				if(relation.equals("=") && !gb.equals("na") && !ref.equals("na")){
					if(!gb.equals(ref)){putAlias(a, gb, e, "GENBANK_EQUAL");}
				}
			}
		}}finally{require(!in.close(), "Assembly report stream read failed");}
	}
	static int reportCompartment(String location, String unit){
		assert(location!=null && unit!=null):"Report columns must be present before biological classification";
		if(location.equalsIgnoreCase("Mitochondrion")){return MITO;}
		if(location.equalsIgnoreCase("Chloroplast") || location.equalsIgnoreCase("Plastid")){return PLASTID;}
		// Explicit Plasmid/other locations are not made nuclear by their assembly unit.
		if(!location.equalsIgnoreCase("Chromosome") && !location.equalsIgnoreCase("Linkage Group")
			&& !location.equalsIgnoreCase("na")){return UNKNOWN;}
		return unit.equals("Primary Assembly")?NUCLEAR:UNKNOWN;
	}
	static void putAlias(Assembly a, String alias, Evidence e, String kind){
		require(!alias.isEmpty() && !alias.equals("NA"), "A report alias must be a real accession");
		require(a.aliases.put(alias, new Match(e, kind))==null, "Ambiguous report accession for "+a.id+": "+alias);
	}

	void loadGbff(String path){
		assert(hashes.contains(path)):"GBFF evidence must retain its original source hash";
		final ByteFile in=RrnaResourceIO.open(path);String id=null, qualifier=null;long length=0;
		boolean recordOpen=false, inSource=false;Assembly selected=null;
		try{for(byte[] line;(line=in.nextLine())!=null;){
			if(EukOrganellePanelInventory.starts(line, "LOCUS ")){
				require(!recordOpen, "GBFF record without terminator: "+path);recordOpen=true;
				id=null;qualifier=null;selected=null;inSource=false;length=Long.parseLong(EukOrganellePanelInventory.word(line, 2));
			}else if(EukOrganellePanelInventory.starts(line, "VERSION ")){
				id=EukOrganellePanelInventory.word(line, 1);selected=byRecord.get(id);
				if(selected!=null && !selected.gbff.equals(path)){selected=null;}
			}else if(EukOrganellePanelInventory.starts(line, "//")){
				require(recordOpen && id!=null && length>0, "GBFF requires a complete versioned record");
				if(selected!=null){
					require(selected.organelle==null, "Duplicate selected GBFF record: "+id);
					selected.organelle=new Evidence(length, organelleCompartment(qualifier), "GBFF_SOURCE", path, hashes.get(path), id, "NA", "NA", qualifier==null?"NA":qualifier);
				}
				recordOpen=false;selected=null;inSource=false;
			}else if(selected!=null){
				if(line.length>5 && line[0]==' ' && line[5]!=' '){
					inSource=EukOrganellePanelInventory.starts(line, "     source          ");
				}else if(inSource && EukOrganellePanelInventory.starts(line, "                     /organelle=")){
					final String raw=EukOrganellePanelInventory.text(line).trim();
					require(raw.startsWith("/organelle=\"") && raw.endsWith("\""), "Expected complete source /organelle qualifier: "+id);
					final String value=raw.substring(12, raw.length()-1);
					require(qualifier==null || qualifier.equals(value), "Conflicting source organelle qualifiers: "+id);qualifier=value;
				}else if(line.length>0 && line[0]!=' '){inSource=false;}
			}
		}require(!recordOpen, "Truncated GBFF source: "+path);
		}finally{require(!in.close(), "Original GBFF read failed");}
	}
	static int organelleCompartment(String value){
		if(value==null){return UNKNOWN;}
		final String lower=value.toLowerCase(Locale.ROOT);
		if(lower.equals("mitochondrion") || lower.startsWith("mitochondrion:")){return MITO;}
		if(lower.equals("chloroplast") || lower.equals("plastid") || lower.startsWith("plastid:")){return PLASTID;}
		return UNKNOWN;
	}
	/** Validate the native tag before using the existing exact accession extractor. */
	static String accession(String seqid, int[] taxid){
		assert(seqid!=null):"Match report aliases against the literal FASTA token, not free-text descriptions";
		taxid[0]=0;
		if(!seqid.startsWith("tid|")){return seqid;}
		final int end=seqid.indexOf('|', 4);
		require(end>4 && end<seqid.length()-1 && seqid.indexOf('|', end+1)<0, "Malformed native tid header: "+seqid);
		long parsed=0;
		for(int i=4;i<end;i++){
			final char c=seqid.charAt(i);require(c>='0' && c<='9', "Malformed native taxid: "+seqid);
			parsed=10*parsed+c-'0';require(parsed<=Integer.MAX_VALUE, "Native taxid overflow: "+seqid);
		}
		require(parsed>0, "Native header taxid must be positive: "+seqid);
		taxid[0]=(int)parsed;
		return Euk18sShredSet.accession(seqid);
	}

	void write(String prefix){
		assert(assemblies.size()>0):"Output covers the complete approved panel, including unknown contigs";
		RrnaResourceIO.fresh(prefix, ".contigs.tsv", ".summary.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".contigs.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder row=new ByteBuilder();long total=0;
		try{
			out.println("selection_group\tassembly\tseqid\tlength\tcompartment\tevidence\tsource\tsource_sha80\tmatched_accession\traw_location\tassembly_unit\tsource_organelle_qualifier\tnative_taxid\tassembly_taxid\ttaxonomy_status");
			summary.println("selection_group\tassembly\tcontigs\tbases\tnuclear_contigs\tmitochondrial_contigs\tplastid_contigs\tunknown_contigs\tnuclear_bases\tmitochondrial_bases\tplastid_bases\tunknown_bases\ttaxonomy_mismatch_contigs");
			for(Assembly a:assemblies){
				if(!a.gbff.equals("NA")){require(a.organelle!=null, "Selected accession absent from original GBFF: "+a.record);}
				final long[] counts=new long[4], bases=new long[4];final ObjectSet<String> seen=new ObjectSet<String>(String.class), seenAccessions=new ObjectSet<String>(String.class);
				final int[] nativeTaxid=new int[1];final long[] mismatches=new long[1];
				RrnaResourceIO.stream(a.fasta, read->{
					final String seqid=RrnaResourceIO.first(read.id);require(seen.add(seqid), "Duplicate FASTA seqid within "+a.id+": "+seqid);
					final String accession=accession(seqid, nativeTaxid);
					// FASTA and assembly metadata can come from different taxonomy snapshots.
					final String taxonomy=nativeTaxid[0]==0?"NOT_TAGGED":a.taxid<=0?"ASSEMBLY_TAXID_UNKNOWN":nativeTaxid[0]==a.taxid?"MATCH":"MISMATCH";
					if(taxonomy.equals("MISMATCH")){mismatches[0]++;}
					require(seenAccessions.add(accession), "Ambiguous normalized FASTA accession within "+a.id+": "+accession);
					final Match match=a.aliases.get(accession);Evidence e=match==null?null:match.evidence;
					String kind=match==null?"UNMAPPED":e.origin+"_"+match.kind;
					if(a.organelle!=null && accession.equals(a.record)){
						require(e==null, "Conflicting report and GBFF evidence sources require adjudication: "+a.id+" "+seqid);
						e=a.organelle;kind=e.compartment==UNKNOWN?"GBFF_SOURCE_ORGANELLE_UNKNOWN":"GBFF_SOURCE_ORGANELLE";
					}
					if(e!=null){require(e.length==read.length(), "FASTA/evidence length mismatch for "+a.id+" "+seqid+": fasta="+read.length()+" evidence="+e.length);}
					final int comp=e==null?UNKNOWN:e.compartment;counts[comp]++;bases[comp]+=read.length();
					row.clear().append(a.group).tab().append(a.id).tab().append(seqid).tab().append(read.length()).tab().append(NAMES[comp]).tab().append(kind);
					if(e==null){row.append("\tNA\tNA\tNA\tNA\tNA\tNA");}
					else{row.tab().append(e.source).tab().append(e.hash).tab().append(e.accession).tab().append(e.location).tab().append(e.unit).tab().append(e.qualifier);}
					out.print(row.tab().append(nativeTaxid[0]).tab().append(a.taxid).tab().append(taxonomy).nl());
				});
				long n=0, b=0;for(int i=0;i<4;i++){n+=counts[i];b+=bases[i];}
				require(n==a.expectedRecords && b==a.expectedBases, "Input inventory changed for "+a.id+": records="+n+"/"+a.expectedRecords+" bases="+b+"/"+a.expectedBases);
				row.clear().append(a.group).tab().append(a.id).tab().append(n).tab().append(b);
				for(long value:counts){row.tab().append(value);}for(long value:bases){row.tab().append(value);}summary.print(row.tab().append(mismatches[0]).nl());total+=n;
			}
		}finally{RrnaResourceIO.finish(out, summary);}
		System.err.println("COMPARTMENTS_COMPLETE assemblies="+assemblies.size()+" contigs="+total+"; no gene truth or annotation-coordinate validation");
	}
	public static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final int NUCLEAR=0, MITO=1, PLASTID=2, UNKNOWN=3;
	static final String[] NAMES={"NUCLEAR", "MITOCHONDRION", "PLASTID", "UNKNOWN"};
	static final class Assembly{
		Assembly(String g, String i, String f, long b, long n, int t){group=g;id=i;fasta=f;expectedBases=b;expectedRecords=n;taxid=t;}
		final String group,id,fasta;final long expectedBases,expectedRecords;final int taxid;String gbff,record;boolean reportSeen;Evidence organelle;
		final ObjectMap<String,Match> aliases=new ObjectMap<String,Match>(String.class, Match.class);
	}
	static final class Evidence{
		Evidence(long l, int c, String o, String s, String h, String a, String loc, String u, String q){length=l;compartment=c;origin=o;source=s;hash=h;accession=a;location=loc;unit=u;qualifier=q;}
		final long length;final int compartment;final String origin,source,hash,accession,location,unit,qualifier;
	}
	static final class Match{Match(Evidence e, String k){evidence=e;kind=k;}final Evidence evidence;final String kind;}
	final ObjectMap<String,String> hashes=new ObjectMap<String,String>(String.class, String.class);
	final ObjectMap<String,Assembly> byAssembly=new ObjectMap<String,Assembly>(String.class, Assembly.class), byRecord=new ObjectMap<String,Assembly>(String.class, Assembly.class);
	final ArrayList<Assembly> assemblies=new ArrayList<Assembly>();
	final ArrayList<String> gbffPaths=new ArrayList<String>();
}
