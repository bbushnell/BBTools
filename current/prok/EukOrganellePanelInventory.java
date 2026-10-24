package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import tax.TaxNode;
import tax.TaxTree;

/** Read-only candidate inventory, never a frozen benchmark or biological truth.
 * Uses existing assembly reports and organelle GBFF metadata; no gene-count selection.
 * @author Ganyu */
public final class EukOrganellePanelInventory {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukOrganellePanelInventory.class,false);
		String cohort=null,reports=null,summary=null,plastid=null,mito=null,gbff=null,tree=null,out=null;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(key.equals("cohort")){cohort=v;}else if(key.equals("reports")){reports=v;}
			else if(key.equals("summary")){summary=v;}else if(key.equals("plastid")){plastid=v;}
			else if(key.equals("mito")){mito=v;}else if(key.equals("gbff")){gbff=v;}
			else if(key.equals("tree")){tree=v;}else if(key.equals("outprefix")){out=v;}
			else{throw new IllegalArgumentException(key);}
		}
		require(cohort!=null && reports!=null && summary!=null && plastid!=null && mito!=null && gbff!=null && tree!=null && out!=null,"All frozen inputs must be named");
		Shared.setThreads(1);
		final EukOrganellePanelInventory x=new EukOrganellePanelInventory();
		x.tree=TaxTree.loadTaxTree(tree,System.err,false,false);
		x.cohort(cohort);x.summary(summary);x.reports(reports);
		x.organelle(plastid,"plastid",15233);x.organelle(mito,"mito",18236);
		for(String file:gbff.split(",")){x.gbff(file);}
		for(Candidate c:x.rows){
			if(!c.comp.equals("nuclear")){require(c.gbffSeen,"Missing GBFF record: "+c.id);}
			x.classify(c);
		}
		x.guardConflicts();x.write(out);Shared.closeStream(pp.outstream);
	}

	void cohort(String path){
		assert(rows.isEmpty()):"One frozen cohort per inventory";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			final String h=text(in.nextLine());require(h.startsWith("index\tclade\tassembly\tbases\tsequence_records\t"),"Unexpected cohort schema");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==20,"Frozen cohort has20 fields");
				final Candidate c=new Candidate("nuclear",p.parseString(2));
				c.sourceClade=p.parseString(1);c.bases=p.parseLong(3);c.records=p.parseLong(4);c.fasta=p.parseString(16);c.hash=p.parseString(17);c.provenance=path;
				if(c.id.startsWith("tid_")){final int end=c.id.indexOf('_',4);require(end>4,"Taxid-prefixed source ID");c.taxid=Integer.parseInt(c.id.substring(4,end));}
				put(c);
			}
		}finally{require(!in.close(),"Cohort input read failed");}
		require(rows.size()==1660,"Preserve every source assembly including zero-rRNA assemblies");
	}

	void summary(String path){
		assert(!rows.isEmpty()):"Only query the declared cohort";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			require(text(in.nextLine()).startsWith("#assembly_accession\ttaxid\tassembly_type\t"),"Unexpected shrunk summary schema");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);final Candidate c=byId.get("nuclear\t"+p.parseString(0));if(c==null){continue;}
				final int extra=p.terms()==12?1:0;require(p.terms()==11+extra,"Unexpected summary columns");
				if(extra==1){require(p.termEquals("",1),"Only the documented empty second field may shift columns");}
				final int tid=p.parseInt(1+extra);c.summaryTaxid=tid;c.summarySource=path;
				if(c.taxid!=0 && c.taxid!=tid){c.taxonomyConflict=true;c.taxid=0;}else{c.taxid=tid;}
			}
		}finally{require(!in.close(),"Summary input read failed");}
	}

	void reports(String path){
		assert(byId.size()==1660):"Assembly reports enrich the frozen cohort; they never select loci";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');Report r=null;
		try{
			for(byte[] line;(line=in.nextLine())!=null;){
				if(starts(line,"# Assembly name:")){apply(r,path);r=new Report();}
				if(r==null){continue;}
				if(starts(line,"# RefSeq assembly accession:")){r.id=value(line);}
				else if(starts(line,"# Taxid:")){r.taxid=Integer.parseInt(value(line));}
				else if(starts(line,"# Assembly level:")){r.level=value(line);}
				else if(starts(line,"# Genome representation:")){r.representation=value(line);}
				else if(line.length>0 && line[0]!='#'){
					p.set(line);require(p.terms()==10,"Assembly report sequence rows require10 columns");
					r.bases+=p.parseLong(8);r.records++;
				}
			}apply(r,path);
		}finally{require(!in.close(),"Assembly report stream read failed");}
	}
	void apply(Report r,String path){
		if(r==null || r.id.equals("NA")){return;}
		final Candidate c=byId.get("nuclear\t"+r.id);if(c==null){return;}
		require(!c.reportSeen,"Duplicate assembly report: "+r.id);c.reportSeen=true;
		c.reportTaxid=r.taxid;
		if(c.taxid!=0 && c.taxid!=r.taxid){c.taxonomyConflict=true;}
		c.taxid=c.taxonomyConflict?0:r.taxid;c.level=r.level;c.representation=r.representation;
		c.reportBases=r.bases;c.reportRecords=r.records;c.qualitySource=path;
	}

	void organelle(String path,String comp,int expected){
		assert(comp.equals("plastid") || comp.equals("mito")):"Organelle identity includes its compartment";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int count=0;
		try{
			RrnaResourceIO.header(in,"accession\ttaxid\torganelle_class\torganism",path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==4,"Organelle index schema");
				final Candidate c=new Candidate(comp,p.parseString(0));c.taxid=p.parseInt(1);c.organelleClass=p.parseString(2);
				c.provenance=path;c.fasta=path.substring(0,path.lastIndexOf('/'))+"/pairs/"+c.id+".fna.gz";put(c);count++;
			}
		}finally{require(!in.close(),"Organelle index input failed");}
		require(count==expected,"Frozen organelle index count changed: "+comp);
	}

	void gbff(String path){
		assert(path!=null):"GBFF metadata source must be explicit";
		final ByteFile in=RrnaResourceIO.open(path);String id=null;long length=0;StringBuilder def=new StringBuilder();boolean inDef=false,segment=false;
		try{
			for(byte[] line;(line=in.nextLine())!=null;){
				if(starts(line,"LOCUS ")){require(id==null,"GBFF record lacks terminator");length=Long.parseLong(word(line,2));def.setLength(0);segment=false;inDef=false;}
				else if(starts(line,"DEFINITION ")){def.append(text(line).substring(12).trim());inDef=true;}
				else if(inDef && starts(line,"            ")){def.append(' ').append(text(line).trim());}
				else{
					inDef=false;
					if(starts(line,"VERSION ")){id=word(line,1);}
					else if(starts(line,"                     /segment=") || starts(line,"                     /chromosome=")){segment=true;}
					else if(starts(line,"//")){
						require(id!=null && length>0,"GBFF record must have accession.version and length");
						Candidate c=byId.get("plastid\t"+id);if(c==null){c=byId.get("mito\t"+id);}
						require(c!=null && !c.gbffSeen,"Missing or duplicate indexed GBFF accession: "+id);
						c.gbffSeen=true;c.bases=length;c.records=1;c.definition=def.toString();c.qualitySource=path;
						final String d=c.definition.toLowerCase(Locale.ROOT);
						c.segment=segment || d.contains("segment") || d.contains("chromosome");
						c.representation=d.contains("complete genome")?"complete_genome_claim":d.contains("complete sequence")?"complete_sequence_claim":"unconfirmed";
						id=null;
					}
				}
			}require(id==null,"Truncated GBFF record");
		}finally{require(!in.close(),"GBFF input failed");}
	}

	void classify(Candidate c){
		assert(c.bases>0 && c.records>0):"Candidate size must come from frozen inventory or GBFF LOCUS";
		// Missing/newer IDs and disagreeing snapshots are retained but never ranked.
		TaxNode n=c.taxid>0?tree.getNode(c.taxid,true):null;boolean fungi=false,animal=false,land=false,green=false,red=false,stram=false,alve=false,euk=false;
		for(int depth=0;n!=null;depth++){
			require(depth<100,"Taxonomy cycle");
			if(n.level==TaxTree.SPECIES){c.species=n.id;}if(n.level==TaxTree.GENUS){c.genus=n.id;}
			fungi|=n.name.equals("Fungi");animal|=n.name.equals("Metazoa");land|=n.name.equals("Embryophyta");
			green|=n.name.equals("Viridiplantae");red|=n.name.equals("Rhodophyta");stram|=n.name.equals("Stramenopiles");alve|=n.name.equals("Alveolata");euk|=n.name.equals("Eukaryota");
			if(n.id==n.pid){break;}n=tree.getNode(n.pid);
		}
		c.clade=fungi?"fungi":animal?"metazoa":land?"landplants":green||red?"green_red_algae":stram||alve?"stramenopiles_alveolates":euk?"other_protists":"UNCLASSIFIED";
		if(c.comp.equals("nuclear")){
			c.bin="nuclear_"+c.clade;
			if(c.reportSeen && c.representation.equals("full") && c.reportBases==c.bases && c.reportRecords==c.records){
				// G11, 2026-10-17: completeness is the gate; size ranks within it.
				c.tier=c.level.equalsIgnoreCase("Complete Genome") || c.level.equals("Chromosome")?0:c.level.equals("Scaffold")?2:c.level.equals("Contig")?3:4;
			}
		}else if(c.comp.equals("plastid")){
			final boolean rare=c.organelleClass.contains("apicoplast") || c.organelleClass.contains("cyanelle");
			c.bin=rare?"plastid_rare":land?"plastid_landplants":green||red?"plastid_green_red_algae":"plastid_other";
		}else{c.bin=animal?"mito_metazoa":fungi?"mito_fungi":green||red?"mito_plants":"mito_protists";}
		// A complete sequence may be one segment or one gene, not the full compartment.
		if(!c.comp.equals("nuclear") && !c.segment){c.tier=c.representation.equals("complete_genome_claim")?0:4;}
		if(!euk || c.species==0 || c.taxonomyConflict){c.tier=4;c.bin=c.comp+"_UNCLASSIFIED";}
		require(!c.taxonomyConflict || c.taxid==0 && c.species==0 && c.bin.endsWith("_UNCLASSIFIED"),
			"Disagreeing taxonomy snapshots must remain explicit and excluded, never silently reconciled");
	}

	/** Conservative breadth guard: both disputed taxonomies need a supported same-genus alternative.
	 * A genus singleton is held for review even when its broad quota bin has other species. */
	void guardConflicts(){
		assert(tree!=null):"Breadth comparisons use the same explicit taxonomy snapshot";
		for(Candidate c:rows){if(c.taxonomyConflict){
			c.summaryGenus=genus(c.summaryTaxid);c.reportGenus=genus(c.reportTaxid);
			for(Candidate other:rows){
				if(other==c || other.taxonomyConflict || other.tier>=4 || !other.comp.equals(c.comp)){continue;}
				if(c.summaryGenus>0 && c.summaryGenus==other.genus){c.summaryAlternatives++;}
				if(c.reportGenus>0 && c.reportGenus==other.genus){c.reportAlternatives++;}
			}
			c.guard=c.summaryAlternatives>0 && c.reportAlternatives>0?"REDUNDANT_DROPPED":"HOLD_COVERAGE_REVIEW";
		}}
	}
	int genus(int tid){
		assert(tree!=null):"Never infer a genus from an organism label";
		TaxNode n=tid>0?tree.getNode(tid,true):null;
		for(int depth=0;n!=null;depth++){
			require(depth<100,"Taxonomy cycle");if(n.level==TaxTree.GENUS){return n.id;}if(n.id==n.pid){break;}n=tree.getNode(n.pid);
		}return 0;
	}
	void write(String prefix){
		assert(rows.size()==1660+15233+18236):"Conserve every source record before ranking";
		RrnaResourceIO.fresh(prefix,".all.tsv",".candidates.tsv",".summary.tsv",".conflicts.tsv",".conflict_clades.tsv");
		writeConflictClades(prefix+".conflict_clades.tsv");
		Collections.sort(rows,ORDER);
		final ByteStreamWriter all=RrnaResourceIO.writer(prefix+".all.tsv"), chosen=RrnaResourceIO.writer(prefix+".candidates.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteStreamWriter conflicts=RrnaResourceIO.writer(prefix+".conflicts.tsv");
		final ByteBuilder b=new ByteBuilder();
		final ObjectSet<String> selectedSpecies=new ObjectSet<String>(String.class);
		try{
			all.println(HEADER);chosen.println(HEADER);summary.println("bin\trequested\tsource_records\tmetadata_supported_records\tselected\tselected_supported\tselected_bases\tlargest_selected\ttop_half_bases\tstatus\tselected_top_completeness_tier\tselected_supported_fragmented\tselected_metadata_only");
			conflicts.println("accession\tsource_clade\tsummary_taxid\treport_taxid\tsummary_genus\treport_genus\tsummary_supported_alternatives\treport_supported_alternatives\trare_guard\treason\treport_stream\tsummary_source");
			for(Candidate c:rows){if(c.taxonomyConflict){
				require(!c.guard.equals("NA"),"Every excluded conflict must have a durable breadth decision");
				b.clear().append(c.id).tab().append(c.sourceClade).tab().append(c.summaryTaxid).tab().append(c.reportTaxid).tab().append(c.summaryGenus).tab().append(c.reportGenus).tab().append(c.summaryAlternatives).tab().append(c.reportAlternatives).tab().append(c.guard).tab().append("TAXONOMY_CONFLICT").tab().append(c.qualitySource).tab().append(c.summarySource).nl();conflicts.print(b);
			}}
			for(Candidate c:rows){row(all,b,c,0);}
			for(int k=0;k<BINS.length;k++){
				final ArrayList<Candidate> candidates=new ArrayList<Candidate>();int supported=0;
				for(Candidate c:rows){if(c.bin.equals(BINS[k])){candidates.add(c);if(c.tier<4){supported++;}}}
				final ObjectSet<String> species=new ObjectSet<String>(String.class),genera=new ObjectSet<String>(String.class);
				final ArrayList<Candidate> selected=new ArrayList<Candidate>();
				// Below-gate tiers fill only an otherwise underfilled quota, for explicit review.
				for(int tier=0;tier<=4;tier++){for(int pass=0;pass<2;pass++){for(Candidate c:candidates){
					if(selected.size()>=QUOTAS[k] || c.tier!=tier || c.species==0 || species.contains(""+c.species) || selectedSpecies.contains(c.comp+'\t'+c.species)){continue;}
					final String genus=c.genus>0?""+c.genus:"species:"+c.species;
					if(pass==0 && genera.contains(genus)){continue;}
					require(!c.taxonomyConflict,"A candidate cannot enter the ranked proposal with unresolved taxonomy");
					species.add(""+c.species);selectedSpecies.add(c.comp+'\t'+c.species);genera.add(genus);selected.add(c);
				}}}
				long bases=0,largest=0,half=0;int good=0,top=0,rank=0;
				for(Candidate c:selected){rank++;row(chosen,b,c,rank);bases+=c.bases;largest=Math.max(largest,c.bases);if(c.tier<4){good++;}if(c.tier==0){top++;}if(rank<=(QUOTAS[k]+1)/2){half+=c.bases;}}
				b.clear().append(BINS[k]).tab().append(QUOTAS[k]).tab().append(candidates.size()).tab().append(supported).tab().append(selected.size()).tab().append(good).tab().append(bases).tab().append(largest).tab().append(half).tab().append("PROPOSAL_ONLY_NOT_FROZEN").tab().append(top).tab().append(good-top).tab().append(selected.size()-good).nl();summary.print(b);
			}
		}finally{RrnaResourceIO.finish(all,chosen,summary,conflicts);}
		System.err.println("INVENTORY_COMPLETE source_records="+rows.size()+" BUSCO=NOT_MEASURED panel_frozen=false");
	}
	/** Both disputed IDs are examined independently; this does not adjudicate either ID. */
	void writeConflictClades(String path){
		assert(tree!=null):"Conflict severity must use the frozen taxonomy, never labels or downloads";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{
			w.println("accession\ttaxid_source\tdisputed_taxid\tpossible_bin\tquota\tclean_species_available\tsupported_clean_species\ttop_complete_clean_species\tseverity");
			for(Candidate c:rows){if(c.taxonomyConflict){
				for(int side=0;side<2;side++){
					final Candidate probe=new Candidate(c.comp,c.id);probe.taxid=side==0?c.summaryTaxid:c.reportTaxid;
					probe.bases=c.bases;probe.records=c.records;classify(probe);
					int quota=-1;for(int k=0;k<BINS.length;k++){if(BINS[k].equals(probe.bin)){quota=QUOTAS[k];}}
					require(quota>=0,"Conflict classification must map to a declared quota bin");
					final ObjectSet<String> clean=new ObjectSet<String>(String.class),supported=new ObjectSet<String>(String.class),complete=new ObjectSet<String>(String.class);
					for(Candidate other:rows){if(!other.taxonomyConflict && other.species>0 && other.bin.equals(probe.bin)){
						clean.add(""+other.species);if(other.tier<4){supported.add(""+other.species);}if(other.tier==0){complete.add(""+other.species);}
					}}
					final String severity=quota==0?"UNRESOLVED_CLADE":clean.size()==0?"EMPTY_CLADE":clean.size()<quota?"UNDERFILLED_CLADE":"GENUS_LOSS_ONLY";
					b.clear().append(c.id).tab().append(side==0?"summary":"report").tab().append(probe.taxid).tab().append(probe.bin).tab().append(quota).tab().append(clean.size()).tab().append(supported.size()).tab().append(complete.size()).tab().append(severity).nl();w.print(b);
				}
			}}
		}finally{RrnaResourceIO.finish(w);}
	}
	void row(ByteStreamWriter w,ByteBuilder b,Candidate c,int rank){
		assert(c.bin!=null):"Classify every row before output";
		b.clear().append(c.bin).tab().append(rank).tab().append(c.comp).tab().append(c.id).tab().append(c.taxid).tab().append(c.species).tab().append(c.genus).tab().append(c.bases).tab().append(c.records).tab().append(c.tier).tab().append(c.level).tab().append(c.representation).tab().append(c.reportBases).tab().append(c.reportRecords).tab().append(c.segment?1:0).tab().append(c.taxonomyConflict?"TAXONOMY_CONFLICT_EXCLUDED":c.tier==0?"TOP_COMPLETENESS_METADATA_SUPPORTED":c.tier<4?"BELOW_GATE_FRAGMENTED_REVIEW":"BELOW_GATE_METADATA_ONLY_REVIEW").tab().append(c.fasta).tab().append(c.hash).tab().append(c.provenance).tab().append(c.qualitySource).tab().append(c.definition).tab().append(c.summaryTaxid).tab().append(c.reportTaxid).nl();w.print(b);
	}
	void put(Candidate c){require(byId.get(c.comp+'\t'+c.id)==null,"Duplicate accession/compartment: "+c.id);byId.put(c.comp+'\t'+c.id,c);rows.add(c);}
	static boolean starts(byte[] b,String s){if(b.length<s.length()){return false;}for(int i=0;i<s.length();i++){if(b[i]!=s.charAt(i)){return false;}}return true;}
	static String text(byte[] b){require(b!=null,"Missing required header");return new String(b,StandardCharsets.UTF_8);}
	static String value(byte[] b){final String s=text(b);require(s.indexOf(':')>=0,"Metadata delimiter");return s.substring(s.indexOf(':')+1).trim();}
	static String word(byte[] b,int target){int n=0;for(int i=0;i<b.length;){while(i<b.length && b[i]<=32){i++;}int start=i;while(i<b.length && b[i]>32){i++;}if(i>start && n++==target){return new String(b,start,i-start,StandardCharsets.UTF_8);}}throw new IllegalArgumentException("Missing header token");}
	static void require(boolean b,String why){if(!b){throw new IllegalArgumentException(why);}}
	static final Comparator<Candidate> ORDER=Comparator.comparingInt((Candidate c)->c.tier).thenComparingLong(c->c.bases).thenComparingLong(c->c.records).thenComparing(c->c.id);
	static final String[] BINS={"nuclear_fungi","nuclear_metazoa","nuclear_landplants","nuclear_green_red_algae","nuclear_stramenopiles_alveolates","nuclear_other_protists","plastid_landplants","plastid_green_red_algae","plastid_other","plastid_rare","mito_metazoa","mito_fungi","mito_plants","mito_protists","nuclear_UNCLASSIFIED","plastid_UNCLASSIFIED","mito_UNCLASSIFIED"};
	static final int[] QUOTAS={30,30,20,15,15,10,16,12,8,4,12,12,8,8,0,0,0};
	static final String HEADER="bin\trank_in_bin\tcompartment\taccession_or_source_id\ttaxid\tspecies_taxid\tgenus_taxid\tbases\trecords\tquality_tier\tassembly_level\trepresentation\treport_bases\treport_records\tsegmented_or_chromosomal\tquality_status\tfasta\tfrozen_fasta_sha80\tidentity_source\tquality_source\tdefinition\tsummary_taxid\treport_taxid";
	static final class Candidate{
		Candidate(String comp_,String id_){comp=comp_;id=id_;}
		final String comp,id;String bin,clade,sourceClade="NA",level="NA",representation="NA",fasta="NA",hash="NA",provenance="NA",qualitySource="NA",definition="NA",organelleClass="NA",summarySource="NA",guard="NA";
		int taxid,summaryTaxid,reportTaxid,species,genus,summaryGenus,reportGenus,summaryAlternatives,reportAlternatives,tier=4;long bases,records,reportBases=-1,reportRecords=-1;boolean reportSeen,gbffSeen,segment,taxonomyConflict;
	}
	static final class Report{String id="NA",level="NA",representation="NA";int taxid;long bases,records;}
	TaxTree tree;
	final ArrayList<Candidate> rows=new ArrayList<Candidate>();
	final ObjectMap<String,Candidate> byId=new ObjectMap<String,Candidate>(String.class,Candidate.class);
}
