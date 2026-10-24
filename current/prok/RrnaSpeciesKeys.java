package prok;

import java.util.ArrayList;
import java.util.Locale;
import java.util.function.IntFunction;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;
import tax.TaxNode;
import tax.TaxTree;

/** Derive construction species keys from saved included truth, never genomes.
 * Original taxids and truth stay immutable. An unresolved assembly is audited
 * and blocks the keyed FASTA, so ambiguity cannot silently exclude its members.
 * @author Raiden
 */
public final class RrnaSpeciesKeys {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaSpeciesKeys.class,false);String truth=null,in=null,treeFile=null,out=null;
		for(String arg:pp.args){final int e=arg.indexOf('=');require(e>0,"Expected flag=value");final String k=arg.substring(0,e).toLowerCase(Locale.ROOT),v=arg.substring(e+1);
			if(k.equals("truth")){truth=v;}else if(k.equals("in")){in=v;}else if(k.equals("tree")){treeFile=v;}else if(k.equals("outprefix")){out=v;}else{throw new IllegalArgumentException(k);}}
		require(truth!=null&&in!=null&&treeFile!=null&&out!=null,"Declare saved truth/construction FASTA, pinned tree and fresh output prefix");
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=false;Shared.setThreads(1);ReadWrite.setZipThreads(1);
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;
		final TaxTree tree=TaxTree.loadTaxTreeFile(treeFile);new RrnaSpeciesKeys().run(truth,in,out,id->tree.getNode(id,true));Shared.closeStream(pp.outstream);
	}
	public void run(String truth,String fasta,String out,IntFunction<TaxNode> lookup){
		RrnaResourceIO.fresh(out,".loci.tsv.gz",".assemblies.tsv",".construction.fa.gz");
		final ArrayList<Locus> loci=new ArrayList<Locus>();final ArrayList<Assembly> assemblies=new ArrayList<Assembly>();
		final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class,Locus.class);
		final ObjectMap<String,Assembly> byAssembly=new ObjectMap<String,Assembly>(String.class,Assembly.class);final IntHashMap speciesCache=new IntHashMap();
		final ByteFile input=RrnaResourceIO.open(truth);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(input,TRUTH_HEADER,truth);for(byte[] row;(row=input.nextLine())!=null;){p.set(row);require(p.terms()==22,"Species keys require all included loci from the22-column truth, before construction filtering");
			final Locus l=new Locus();l.id=p.parseString(0);l.clade=p.parseString(1);l.assembly=p.parseString(2);l.taxid=p.parseInt(3);l.seq=p.parseString(4);l.start=p.parseInt(5);l.stop=p.parseInt(6);l.strand=p.parseString(7);l.length=p.parseInt(8);l.comp=p.parseString(20);
			require(l.taxid>0&&l.start>0&&l.stop>=l.start&&l.length==l.stop-l.start+1&&(l.strand.equals("+")||l.strand.equals("-")),"Truth coordinates and original taxonomy must remain valid");
			Euk18sM1Grade.validComp(l.comp);require((p.termEquals("0",17)||p.termEquals("1",17))&&(p.termEquals("0",21)||p.termEquals("1",21)),"Admission and construction flags must be explicit");
			l.construction=p.termEquals("1",21);require(l.construction==(p.termEquals("1",17)&&l.comp.equals("no_keyword")),"Construction flags must match the declared admitted/no_keyword recipe");
			require(l.id.startsWith("tid|"+l.taxid+'|')&&byId.get(l.id)==null,"Original locus IDs must be unique and retain their source taxid");
			int species=speciesCache.get(l.taxid);if(species<0){species=species(l.taxid,lookup);speciesCache.put(l.taxid,species);}l.species=species;
			final String key=l.clade+'\t'+l.assembly;Assembly a=byAssembly.get(key);if(a==null){a=new Assembly(l.clade,l.assembly);byAssembly.put(key,a);assemblies.add(a);}l.group=a;
			a.total++;if(species==0){a.missing++;}else{a.counts.increment(species);if(l.comp.equals("no_keyword")){a.nuclear.increment(species);}}if(l.construction){a.construction++;}
			byId.put(l.id,l);loci.add(l);
		}}finally{require(!input.close(),"Saved truth read failed");}require(!loci.isEmpty(),"An empty truth pool cannot produce a species selection population");
		int unresolved=0;for(Assembly a:assemblies){choose(a);if(a.key==0){unresolved++;}}
		// Write complete ambiguity evidence before refusing to manufacture a key.
		final ByteStreamWriter audit=RrnaResourceIO.writer(out+".loci.tsv.gz"),groups=RrnaResourceIO.writer(out+".assemblies.tsv");final ByteBuilder b=new ByteBuilder();
		try{audit.println("locus\tclade\tassembly\tsource_taxid\tsource_species_taxid\tassembly_species_key\tcompartment\tconstruction_pool\tstatus");
			for(Locus l:loci){audit.print(b.clear().append(l.id).tab().append(l.clade).tab().append(l.assembly).tab().append(l.taxid).tab().append(l.species).tab().append(l.group.key).tab().append(l.comp).tab().append(l.construction?1:0).tab().append(l.group.status).nl());}
			groups.println("clade\tassembly\tincluded\tconstruction\tmissing_species_ancestry\tdistinct_species\tassembly_species_key\tstatus");
			for(Assembly a:assemblies){groups.print(b.clear().append(a.clade).tab().append(a.id).tab().append(a.total).tab().append(a.construction).tab().append(a.missing).tab().append(a.counts.size()).tab().append(a.key).tab().append(a.status).nl());}
		}finally{RrnaResourceIO.finish(audit,groups);}
		require(unresolved==0,"Unresolved species ancestry or tie in "+unresolved+" assemblies; inspect saved audits. No keyed FASTA emitted.");
		final ByteStreamWriter output=RrnaResourceIO.writer(out+".construction.fa.gz");final long[] count={0};
		try{RrnaResourceIO.stream(fasta,r->{final Locus l=byId.get(RrnaResourceIO.first(r.id));require(l!=null&&l.construction&&!l.seen,"Construction FASTA must contain each declared construction member exactly once");
			for(String field:new String[]{"source_taxid","source_species_taxid","assembly_species_key","source_locus"}){require(RrnaResourceIO.optionalField(r.id,field)==null,"Input must be original construction FASTA, not already keyed");}
			require(r.length()==l.length&&RrnaResourceIO.field(r.id,"assembly").equals(l.assembly)&&RrnaResourceIO.field(r.id,"seqid").equals(l.seq)&&Integer.parseInt(RrnaResourceIO.field(r.id,"start"))==l.start&&Integer.parseInt(RrnaResourceIO.field(r.id,"stop"))==l.stop&&RrnaResourceIO.field(r.id,"strand").equals(l.strand)&&RrnaResourceIO.field(r.id,"compartment").equals(l.comp),"Saved FASTA must match the exact truth locus and metadata");
			final int delimiter=r.id.indexOf('|',4);require(delimiter>4,"Original tid header required");
			output.print(b.clear().append(">tid|").append(l.group.key).append('|').append(r.id.substring(delimiter+1)).append(" source_taxid=").append(l.taxid).append(" source_species_taxid=").append(l.species).append(" assembly_species_key=").append(l.group.key).append(" source_locus=").append(l.id).nl().append(r.bases).nl());l.seen=true;count[0]++;
		});}finally{RrnaResourceIO.finish(output);}
		long expected=0;for(Locus l:loci){if(l.construction){require(l.seen,"Missing construction member "+l.id);expected++;}}
		require(count[0]==expected,"Species rekeying must conserve every construction member without selection");
		System.err.println("RRNA_SPECIES_KEYS_PASS included="+loci.size()+" assemblies="+assemblies.size()+" construction="+expected+" original_truth_unchanged no_genome_scan");
	}
	public static int species(int source,IntFunction<TaxNode> lookup){
		require(source>0&&lookup!=null,"Resolve a declared source TaxID in the pinned taxonomy");final IntHashMap seen=new IntHashMap(16);TaxNode n=lookup.apply(source);
		while(n!=null){require(seen.get(n.id)<0,"Cycle in taxonomy ancestry for "+source);seen.put(n.id,1);
			// Original rank, not a rank inferred by TaxTree simplification, is NCBI's rank.
			if(n.originalLevel()==TaxTree.SPECIES_E){return n.id;}if(n.pid==n.id){return 0;}n=lookup.apply(n.pid);
		}return 0;
	}
	public static void choose(Assembly a){
		require(a.total>0,"Assembly key is a vote over included loci, never a filename guess");if(a.missing>0){a.status="MISSING_SPECIES_ANCESTRY";return;}
		int max=-1,bestNuclear=-1,winner=0,ties=0,primaryTies=0;
		for(int key:a.counts.keys()){if(key==a.counts.invalid()){continue;}final int n=a.counts.get(key),nu=Math.max(0,a.nuclear.get(key));
			if(n>max){max=n;bestNuclear=nu;winner=key;ties=1;primaryTies=1;}else if(n==max){primaryTies++;if(nu>bestNuclear){bestNuclear=nu;winner=key;ties=1;}else if(nu==bestNuclear){ties++;}}
		}
		if(ties!=1){a.status="UNRESOLVED_TIE";return;}a.key=winner;a.status=a.counts.size()==1?"ONE_SPECIES":primaryTies>1?"TIE_NO_KEYWORD":"MAJORITY_INCLUDED";
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Locus{String id,clade,assembly,seq,strand,comp;int taxid,species,start,stop,length;boolean construction,seen;Assembly group;}
	public static final class Assembly{public Assembly(String c,String a){clade=c;id=a;}final String clade,id;public final IntHashMap counts=new IntHashMap(4),nuclear=new IntHashMap(4);public int total,construction,missing,key;public String status;}
	public static final String TRUTH_HEADER="locus\tclade\tassembly\ttaxid\tseqid\tstart\tstop\tstrand\tlength\tcontig_length\tcontig_edge\tundefined\tmodel_from\tmodel_to\ttrunc\tcm_score\tevalue\tfull_length_admitted\torganelle_keyword\ttable_line\tcompartment\tconstruction_pool";
}
