package prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Reports original GFF context at exact prediction coordinates; never regrades or admits truth.
 * The execution producer must verify both supplied sha80 manifests before invocation.
 * All source rows that overlap are retained, including parents and opposite strands.
 * Circular features use explicit region topology; their original coordinates remain unchanged.
 * @author Brian, Ganyu
 */
public final class NcrnaReviewFeatureContext {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaReviewFeatureContext.class, false);
		String loci=null, roster=null, gffHashes=null, fastaHashes=null, out=null, family=null, model=null;
		int expected=-1, genomes=-1, contigs=-1, modelLength=-1;
		try{
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');require(pair.length==2, "Use flag=value");
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair[1];
				if(key.equals("loci")){loci=value;}else if(key.equals("roster")){roster=value;}
				else if(key.equals("gffhashes")){gffHashes=value;}else if(key.equals("fastahashes")){fastaHashes=value;}
				else if(key.equals("outprefix")){out=value;}else if(key.equals("family")){family=value;}
				else if(key.equals("modelid")){model=value;}else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}
				else if(key.equals("genomes")){genomes=Parse.parseIntKMG(value);}else if(key.equals("contigs")){contigs=Parse.parseIntKMG(value);}
				else if(key.equals("modellength")){modelLength=Parse.parseIntKMG(value);}else{throw new IllegalArgumentException(key);}
			}
			require(loci!=null && roster!=null && gffHashes!=null && fastaHashes!=null && out!=null
				&& family!=null && model!=null && model.matches("RF[0-9]{5}") && expected>0 && genomes>0 && contigs>0 && modelLength>0,
				"Bind prediction metadata, original roster, externally verified GFF/FASTA hash manifests, family/model and denominators");
			Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
			final NcrnaReviewFeatureContext x=new NcrnaReviewFeatureContext(family, model, modelLength);
			x.loadRoster(roster, gffHashes, fastaHashes, genomes);x.loadLoci(loci, expected);
			x.loadContigs(contigs);x.readAnnotations();x.write(out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaReviewFeatureContext(String family_, String model_, int length_){family=family_;model=model_;modelLength=length_;}

	/*--------------------------------------------------------------*/
	/*----------------         Source binding       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadRoster(String path, String gffHashes, String fastaHashes, int expected){
		final ObjectMap<String,String> gh=EukPanelGffValidation.hashes(gffHashes), fh=EukPanelGffValidation.hashes(fastaHashes);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> gffs=new ObjectSet<String>(String.class), fastas=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, NcrnaCalibrationPrep.ROSTER_HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==11, "Original genome roster schema changed");
				final Genome g=new Genome(p.parseString(1), p.parseString(2), p.parseString(6), p.parseString(7));
				g.gffHash=gh.get(g.gff);g.fastaHash=fh.get(g.fasta);
				require(!g.id.isEmpty() && !g.assembly.isEmpty() && g.gffHash!=null && g.fastaHash!=null
					&& byGenome.put(g.id, g)==null && gffs.add(g.gff) && fastas.add(g.fasta),
					"Unique genome/GFF/FASTA bindings with verified source hashes required: "+g.id);
				genomes.add(g);
			}
		}finally{require(!in.close(), "Source roster read failed");}
		require(genomes.size()==expected, "Original genome denominator changed");
	}
	private void loadLoci(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, NcrnaReviewUnion.LOCUS_HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==14 && p.termEquals(family, 5) && p.termEquals("PREDICTION_UNRESOLVED", 13),
					"Consume the explicit prediction union, never an inferred annotation schema");
				final Genome g=byGenome.get(p.parseString(1));
				require(g!=null && p.termEquals(g.assembly, 2) && p.termEquals(g.fasta, 3), "Foreign prediction source binding");
				final Locus l=new Locus(p, g);
				require(ids.add(l.id), "Repeated review ID: "+l.id);loci.add(l);
				Group group=wanted.get(l.seqid);
				if(group==null){group=new Group(g, l.length);wanted.put(l.seqid, group);}
				require(group.genome==g && group.length==l.length, "Ambiguous prediction contig ownership or length");group.loci.add(l);
			}
		}finally{require(!in.close(), "Prediction metadata read failed");}
		require(loci.size()==expected, "Prediction union denominator changed");
	}
	private void loadContigs(int expected){
		assert(contigs.isEmpty()):"Read the original FASTA population exactly once before annotation binding";
		for(Genome g:genomes){
			RrnaResourceIO.stream(g.fasta, r->{
				final String name=RrnaResourceIO.first(r.id);
				require(r.bases!=null && r.length()>0 && !contigs.contains(name), "Missing sequence or ambiguous source contig: "+name);
				final Contig c=new Contig(g, r.length());contigs.put(name, c);
				final Group group=wanted.get(name);
				if(group!=null){
					require(group.genome==g && group.length==r.length() && group.contig==null, "Prediction contig must bind its original FASTA and actual length: "+name);
					group.contig=c;c.group=group;
				}
			});
		}
		require(contigs.size()==expected, "Original contig denominator changed");
		for(Locus l:loci){require(wanted.get(l.seqid).contig!=null, "Prediction contig missing from original FASTAs: "+l.seqid);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Annotation context     ----------------*/
	/*--------------------------------------------------------------*/

	private void readAnnotations(){
		assert(!genomes.isEmpty() && !contigs.isEmpty()):"Original FASTA ownership is established before any GFF row is interpreted";
		for(Genome g:genomes){
			final ByteFile in=RrnaResourceIO.open(g.gff);final LineParser1 p=new LineParser1('\t');long line=0;
			String lastSeq=null;Contig current=null;
			try{
				for(byte[] row;(row=in.nextLine())!=null;){
					line++;if(row.length==0){continue;}
					if(row[0]=='#'){if(Tools.startsWith(row, "##FASTA")){break;}continue;}
					p.set(row);if(p.terms()!=9){throw new IllegalArgumentException("Original GFF9 required: "+g.gff+":"+line);}
					if(lastSeq==null || !p.termEquals(lastSeq, 0)){lastSeq=p.parseString(0);current=contigs.get(lastSeq);}
					final Contig c=current;
					if(c==null || c.genome!=g){throw new IllegalArgumentException("GFF contig absent or foreign to its bound source: "+g.gff+":"+line);}
					final int start=p.parseInt(3), stop=p.parseInt(4);
					if(start<1 || start>c.length || stop<start || (long)stop-start+1>c.length){
						throw new IllegalArgumentException("Source feature must span at most one physical contig turn: "+g.gff+":"+line);
					}
					final boolean region=p.termEquals("region", 2);
					if(region){
						final String circular=RrnaResourceIO.optionalField(p.parseString(8), "Is_circular");
						if(circular!=null){
							require(start==1 && stop==c.length && (circular.equals("true") || circular.equals("false")), "Topology must describe the full actual contig");
							final boolean value=circular.equals("true");require(!c.topologyKnown || c.circular==value, "Conflicting circular topology");c.topologyKnown=true;c.circular=value;
						}
					}
					if(stop>c.length){c.wrapped=true;}
					if(c.group==null){continue;}
					Feature feature=null;
					for(Locus l:c.group.loci){
						final int overlap=overlap(l.start, l.stop, start, stop, c.length);
						if(overlap==0){continue;}
						if(feature==null){feature=new Feature(row, p, line, model);}
						l.overlaps.add(new Overlap(feature, overlap));
					}
				}
			}finally{require(!in.close(), "Original GFF read failed: "+g.gff);}
		}
		for(Contig c:contigs.values()){if(c!=null){require(!c.wrapped || c.topologyKnown && c.circular,
			"Origin-spanning annotation requires explicit circular source topology");}}
	}
	/** Physical overlap of a linear prediction with one possibly origin-spanning source row. */
	public static int overlap(int start, int stop, int sourceStart, int sourceStop, int length){
		assert(start>0 && stop>=start && stop<=length && sourceStart>0 && sourceStart<=length
			&& sourceStop>=sourceStart && (long)sourceStop-sourceStart+1<=length):"Intervals follow PanelAnnotationGrade's single-turn source convention";
		int n=Math.max(0, Math.min(stop, Math.min(sourceStop, length))-Math.max(start, sourceStart)+1);
		if(sourceStop>length){n+=Math.max(0, Math.min(stop, sourceStop-length)-start+1);}return n;
	}
	private static boolean hasRfam(String attrs, String model){
		assert(attrs!=null && model!=null):"Target annotation uses an exact Rfam Dbxref token, not a product substring";
		final String value=RrnaResourceIO.optionalField(attrs, "Dbxref");if(value==null){return false;}
		for(int from=0;from<value.length();){
			int to=value.indexOf(',', from);if(to<0){to=value.length();}
			if(to-from==5+model.length() && value.regionMatches(true, from, "Rfam:", 0, 5)
				&& value.regionMatches(from+5, model, 0, model.length())){return true;}from=to+1;
		}return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String prefix){
		RrnaResourceIO.fresh(prefix, ".overlaps.tsv", ".loci.tsv", ".summary.tsv");
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+".overlaps.tsv"), rows=RrnaResourceIO.writer(prefix+".loci.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();long total=0, regions=0;int intergenic=0, parentOnly=0;
		try{
			detail.println("review_id\tgenome\tassembly\tfasta\tfasta_sha80\tgff\tgff_sha80\tgff_line\tfeature_class\ttarget_rfam_tag\toverlap_bp\tstrand_relation\torigin_spanning\tseqid\tsource\ttype\tstart1\tstop1\tscore\tstrand\tphase\toriginal_attributes");
			rows.println("review_id\tgenome\tassembly\tseqid\tstart1\tstop1\tstrand\tcall_length\tmodel_length\tlength_minus_model\tlength_ratio\tfeature_rows\tregion_rows\tparent_rows\ttarget_annotation_rows\tother_ncrna_rows\tcds_rows\tother_rows\tsame_strand_rows\topposite_strand_rows\tunknown_strand_rows\tcontext\tgff\tgff_sha80\tfasta\tfasta_sha80");
			for(Locus l:loci){
				final Genome g=l.genome;int features=0, regionCount=0, parents=0, target=0, otherRna=0, cds=0, other=0, same=0, opposite=0, unknown=0;
				for(Overlap o:l.overlaps){
					final Feature f=o.feature;final String relation=f.strand=='.'||f.strand=='?'?"UNKNOWN":f.strand==l.strand?"SAME":"OPPOSITE";
					if(f.kind.equals("REGION_METADATA")){regionCount++;regions++;}
					else{
						features++;total++;
						if(f.kind.equals("PARENT_ANNOTATION")){parents++;}else if(f.kind.equals("TARGET_RNA_ANNOTATION")){target++;}
						else if(f.kind.equals("OTHER_NCRNA")){otherRna++;}else if(f.kind.equals("CDS")){cds++;}else{other++;}
						if(relation.equals("SAME")){same++;}else if(relation.equals("OPPOSITE")){opposite++;}else{unknown++;}
					}
					detail.print(b.clear().append(l.id).tab().append(g.id).tab().append(g.assembly).tab().append(g.fasta).tab().append(g.fastaHash)
						.tab().append(g.gff).tab().append(g.gffHash).tab().append(f.line).tab().append(f.kind).tab().append(f.target).tab().append(o.bases)
						.tab().append(relation).tab().append(f.stop>l.length).tab().append(f.raw).nl());
				}
				require(features==parents+target+otherRna+cds+other && features==same+opposite+unknown, "Every overlapping source row must retain type and strand accounting");
				final String context=features==0?"INTERGENIC":features==parents?"PARENT_ONLY_OR_AMBIGUOUS":"FEATURE_OVERLAP";
				if(features==0){intergenic++;}else if(features==parents){parentOnly++;}
				final int length=l.stop-l.start+1;
				rows.print(b.clear().append(l.id).tab().append(g.id).tab().append(g.assembly).tab().append(l.seqid).tab().append(l.start).tab().append(l.stop)
					.tab().append(l.strand).tab().append(length).tab().append(modelLength).tab().append(length-modelLength).tab().append((double)length/modelLength, 6)
					.tab().append(features).tab().append(regionCount).tab().append(parents).tab().append(target).tab().append(otherRna).tab().append(cds).tab().append(other)
					.tab().append(same).tab().append(opposite).tab().append(unknown).tab().append(context).tab().append(g.gff).tab().append(g.gffHash).tab().append(g.fasta).tab().append(g.fastaHash).nl());
			}
			summary.println("loci\tgenomes\tcontigs\tfeature_overlap_rows\tregion_overlap_rows\tintergenic_loci\tparent_only_loci\tinterpretation");
			summary.print(b.clear().append(loci.size()).tab().append(genomes.size()).tab().append(contigs.size()).tab().append(total).tab().append(regions)
				.tab().append(intergenic).tab().append(parentOnly).append("\tSOURCE_CONTEXT_NOT_CM_ADMISSION_OR_FP_VERDICT\n"));
		}finally{RrnaResourceIO.finish(detail, rows, summary);}
		System.err.println("NCRNA_REVIEW_FEATURE_CONTEXT_PASS loci="+loci.size()+" overlapRows="+total+" intergenic="+intergenic+" parentOnly="+parentOnly);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome {
		Genome(String id_, String assembly_, String gff_, String fasta_){id=id_;assembly=assembly_;gff=gff_;fasta=fasta_;}
		final String id, assembly, gff, fasta;String gffHash, fastaHash;
	}
	private static final class Contig {
		Contig(Genome genome_, int length_){genome=genome_;length=length_;}
		final Genome genome;final int length;boolean topologyKnown, circular, wrapped;Group group;
	}
	private static final class Group {
		Group(Genome genome_, int length_){genome=genome_;length=length_;}
		final Genome genome;final int length;Contig contig;final ArrayList<Locus> loci=new ArrayList<Locus>();
	}
	private static final class Locus {
		Locus(LineParser1 p, Genome g){
			id=p.parseString(0);seqid=p.parseString(4);genome=g;final String s=p.parseString(6);
			require(!id.isEmpty() && !seqid.isEmpty() && (s.equals("+") || s.equals("-")), "Explicit prediction ID/contig/strand required");strand=s.charAt(0);
			start=p.parseInt(7);stop=p.parseInt(8);length=p.parseInt(9);final int left=p.parseInt(10), right=p.parseInt(11), context=p.parseInt(12);
			require(start>0 && stop>=start && stop<=length && left>=0 && right>=0 && (long)context==(long)stop-start+1+left+right,
				"Prediction coordinates and actual context dimensions must agree");
			require(strand=='+'?start-left>=1 && (long)stop+right<=length:start-right>=1 && (long)stop+left<=length, "Prediction context must remain inside its bound contig");
		}
		final String id, seqid;final Genome genome;final char strand;final int start, stop, length;
		final ArrayList<Overlap> overlaps=new ArrayList<Overlap>();
	}
	private static final class Feature {
		Feature(byte[] row, LineParser1 p, long line_, String model){
			line=line_;stop=p.parseInt(4);final String s=p.parseString(6), type=p.parseString(2);
			require(s.equals("+") || s.equals("-") || s.equals(".") || s.equals("?"), "GFF strand must be explicit or unknown");strand=s.charAt(0);
			final String attrs=p.parseString(8);target=hasRfam(attrs, model);
			kind=type.equals("region")?"REGION_METADATA":type.equals("gene")||type.equals("pseudogene")||type.equals("mRNA")||type.equals("transcript")?"PARENT_ANNOTATION"
				:type.equals("CDS")?"CDS":target?"TARGET_RNA_ANNOTATION":type.endsWith("RNA")||type.endsWith("_RNA")||RrnaResourceIO.optionalField(attrs, "ncRNA_class")!=null?"OTHER_NCRNA":"OTHER_FEATURE";
			raw=row.clone(); // ByteFile implementations may reuse their input buffer; retained evidence owns its bytes.
		}
		final byte[] raw;final String kind;final long line;final int stop;final char strand;final boolean target;
	}
	private static final class Overlap {
		Overlap(Feature feature_, int bases_){feature=feature_;bases=bases_;}
		final Feature feature;final int bases;
	}
	private final String family, model;private final int modelLength;
	private final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ObjectMap<String,Group> wanted=new ObjectMap<String,Group>(String.class, Group.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();private final ArrayList<Locus> loci=new ArrayList<Locus>();
}
