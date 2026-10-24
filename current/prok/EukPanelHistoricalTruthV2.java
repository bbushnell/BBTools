package prok;

import java.util.ArrayList;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import structures.ByteBuilder;
import prok.EukPanelHistoricalDisposition.Group;
import prok.EukPanelHistoricalDisposition.Locus;
import prok.EukPanelHistoricalDisposition.Row;

/** Materializes the reviewed class rule without freezing or rewriting truth_v1.
 * Existing same-class loci retain every biological field. Unresolved duplicate,
 * promotion or new ribosomal crossfold cases stop before materialization.
 * @author Brian, Raiden */
final class EukPanelHistoricalTruthV2 {
	EukPanelHistoricalTruthV2(EukPanelHistoricalDisposition source,String inputs,int expected){x=source;roster=inputs;expectedAssemblies=expected;}
	void write(String history, String prefix){
		assert(!x.rows.isEmpty() && !x.baselineRows.isEmpty()):"Materialization follows the complete evidence and immutable-baseline join";
		prepare();
		RrnaResourceIO.fresh(prefix, ".truth_candidate.tsv.gz", ".historical_masks.tsv.gz", ".materialization.tsv");
		final ByteStreamWriter truth=RrnaResourceIO.writer(prefix+".truth_candidate.tsv.gz"), masks=RrnaResourceIO.writer(prefix+".historical_masks.tsv.gz"), counts=RrnaResourceIO.writer(prefix+".materialization.tsv");
		final ByteBuilder b=new ByteBuilder();final LineParser1 p=new LineParser1('\t');
		try{
			truth.println(EukPanelTruthFreeze.HEADER);
			for(Locus l:x.baselineRows){
				p.set(l.raw);b.clear();
				for(int i=0; i<35; i++){if(i>0){b.tab();}
					if(i==32){b.append(EukPanelTruthFreeze.V2_CANDIDATE);}else if(i==34){b.append(EukPanelTruthFreeze.V2_POLICY);}else{field(b,p,i);}
				}
				truth.print(b.nl());increment(l.family,0);
			}
			for(Row r:additions){newLocus(truth,b,p,r);increment(r.family,r.hit.ga?1:2);}
			masks.println(EukPanelHistoricalMask.OUTPUT_HEADER);final long[] seen={0,0};
			EukPanelHistoricalEvidence.read(history,EukPanelHistoricalMask.OUTPUT_HEADER,q->{
				seen[0]++;final Row r=x.bySource.get(q.parseString(0)+'\t'+q.parseString(8)+'\t'+q.parseLong(10));
				if(r==null){masks.println(q.line());return;}seen[1]++;
				final boolean retained=r.action.equals("BELOW_GA");
				final String category=retained?"HIST_LINKED_UNSCORED":r.action.equals("ADMITTED")?"ADMITTED_LINK":r.opposite!=null?"ORIENTATION_ERROR":"UNSUPPORTED_REVIEWED";
				b.clear();for(int i=0; i<14; i++){if(i>0){b.tab();}
					if(i==6){b.append(category);}else if(i==7){b.append(retained?1:0);}
					else if(i==13){b.append("HISTORICAL_CLASS_REVIEW:").append(r.id).append(':').append(r.action);}
					else{field(b,q,i);}
				}masks.print(b.nl());increment(r.family,retained?4:5);
			});
			require(seen[1]==x.rows.size(),"Every reviewed historical mask must receive exactly one disposition");
			counts.println("fold_family\tbaseline_loci\tadded_admitted\tadded_below_ga\tlinked_existing\tretained_historical_masks\tlifted_historical_masks\tcombined_loci\tstatus");
			long total=0;for(String family:EukPanelCmProposal.sorted(totals)){
				final long[] n=totals.get(family);b.clear().append(family);for(long v:n){b.tab().append(v);}
				final long combined=n[0]+n[1]+n[2];total+=combined;
				counts.print(b.tab().append(combined).tab().append(EukPanelTruthFreeze.V2_CANDIDATE).nl());
			}
			require(total==x.baselineRows.size()+additions.size(),"Materialized truth must conserve baseline plus unique additions");
			System.err.println("HISTORICAL_TRUTH_V2_CANDIDATE baseline="+x.baselineRows.size()+" added="+additions.size()+" total="+total+" history="+seen[0]+" frozen=false");
		}finally{RrnaResourceIO.finish(truth,masks,counts);}
	}
	void prepare(){
		assert(additions.isEmpty()):"Resolve every physical-locus decision before opening candidate output files";
		EukPanelHistoricalEvidence.read(roster,INPUT_HEADER,p->{
			require(p.terms()==15 && (p.termEquals("nuclear",0) || p.termEquals("mito",0) || p.termEquals("plastid",0)),"Complete panel roster schema/compartment changed");
			final String assembly=p.parseString(1),bin=p.parseString(2);
			require(!assembly.isEmpty() && !bin.isEmpty() && !bin.equals("NA") && bins.put(assembly,bin)==null,"Panel assembly/bin must be unique and explicit");
		});
		require(bins.size()==expectedAssemblies,"Panel roster must retain assemblies with zero baseline loci");
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		for(Locus l:x.baselineRows){
			require(ids.add(l.id),"Duplicate baseline ID");
			require(l.bin.equals(bins.get(l.assembly)),"Baseline assembly selection bin disagrees with complete panel roster");
		}
		for(Group g:x.groups){
			require(g.rows.size()==1 && g.overlapping.isEmpty(),"Duplicate or overlapping proposed groups require explicit locus review: "+g.id);
			final Row r=g.seed;
			if(!r.sameLinks.isEmpty()){
				require(r.sameLinks.size()==1 && r.sameLinks.get(0).admitted==r.hit.ga,"Multiple baseline links or a changed admission require explicit review: "+r.id);
				increment(r.family,3);continue;
			}
			require(ids.add(g.id) && bins.get(r.assembly)!=null,"New locus requires a unique ID and an established assembly selection bin: "+r.id);
			if(r.hit.ga && ribosomal(r.family)){
				for(Locus l:x.baselineRows){
					if(l.admitted && ribosomal(l.family) && !l.family.equals(r.family) && l.assembly.equals(r.assembly) && l.seqid.equals(r.seqid) && l.strand.equals(r.strand)){
						require(overlap(l.start,l.end,r.hit.start,r.hit.end)==0,"New ribosomal crossfold overlap needs explicit boundary-eligibility review: "+r.id+'/'+l.id);
					}
				}
				for(Row other:additions){
					if(other.hit.ga && ribosomal(other.family) && !other.family.equals(r.family) && other.assembly.equals(r.assembly) && other.seqid.equals(r.seqid) && other.strand.equals(r.strand)){
						require(overlap(other.hit.start,other.hit.end,r.hit.start,r.hit.end)==0,"New crossfold additions need explicit review");
					}
				}
			}
			additions.add(r);
		}
	}
	void newLocus(ByteStreamWriter out,ByteBuilder b,LineParser1 p,Row r){
		assert(r.sameLinks.isEmpty()):"Only unmatched supported evidence may create a new reference locus";
		p.set(r.raw);final boolean ga=r.hit.ga,full=p.parseInt(84)==1 && p.parseInt(85)==p.parseInt(86),truncated=!p.termEquals("no",87);
		final String edge=r.hit.start==1?(r.hit.end==r.contigLength?"BOTH":"LEFT"):r.hit.end==r.contigLength?"RIGHT":"INTERNAL";
		final String geometry=truncated || !edge.equals("INTERNAL")?"EDGE":full?"FULLSPAN_INTACT":"PARTIAL_INTERNAL";
		final String model=p.parseString(79),flags="HIST_SUPPORTED;HIST_ANY_RNA;HIST_SAME_FAMILY";
		b.clear().append(r.group.id).tab().append(r.assembly).tab().append(r.seqid).tab().append(r.strand).tab().append(r.hit.start).tab().append(r.hit.end)
			.tab().append(r.family).tab().append(bins.get(r.assembly)).tab();field(b,p,16);
		b.tab().append(model).tab().append("HISTORICAL_").append(r.id).tab();field(b,p,28);
		b.tab();field(b,p,80);b.tab();field(b,p,88);
		b.tab().append(ga?1:0).tab().append(full?1:0).tab().append(truncated?1:0).tab().append(edge).append("\t1\t").append(model)
			.append("\t0\t0\tNONE\tREVIEWED_CLASS_ADDENDUM\t").append(ga?1:0).tab().append(ga?model:"NA").tab().append(ga?"ADMITTED":"BELOW_GA")
			.tab().append(geometry).tab().append(flags).append("\t1\t1\t").append(ga?1:0).tab().append(EukPanelTruthFreeze.V2_CANDIDATE)
			.tab().append(ga?1:0).tab().append(EukPanelTruthFreeze.V2_POLICY);
		out.print(b.nl());
	}
	void increment(String family,int column){
		assert(column>=0 && column<6):"Materialization counters distinguish baseline, additions, links and mask disposition";
		long[] n=totals.get(family);if(n==null){n=new long[6];totals.put(family,n);}n[column]++;
	}
	static void field(ByteBuilder b,LineParser1 p,int i){p.setBounds(i);b.append(p.line(),p.a(),p.b()-p.a());}
	static boolean ribosomal(String family){return family.equals("SSU") || family.equals("LSU") || family.equals("5S") || family.equals("5.8S");}
	static long overlap(long a,long b,long c,long d){return Math.max(0,Math.min(b,d)-Math.max(a,c)+1);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	final EukPanelHistoricalDisposition x;
	final String roster;final int expectedAssemblies;
	static final String INPUT_HEADER="compartment\tassembly\tbin\ttaxid\tbases\tsequence_records\tfasta\tselected_fasta_sha80\tgff\tannotation_status\tgbff\tgbff_record_accession\tmetadata_status\tannotation_source_fasta\tpairing_provenance";
	final ObjectMap<String,String> bins=new ObjectMap<String,String>(String.class,String.class);
	final ObjectMap<String,long[]> totals=new ObjectMap<String,long[]>(String.class,long[].class);
	final ArrayList<Row> additions=new ArrayList<Row>();
}
