package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Original annotation overlaps and nearest downstream CDS for frozen B12 calls.
 * The producer verifies source GFF sha80s before invocation. Gene/product hints
 * are descriptive, never CM support, truth admission or a biological verdict.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminGenomicContext {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminGenomicContext.class, false);
		String loci=null, context=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("loci")){loci=v;}else if(k.equals("context")){context=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}
		}
		require(loci!=null && context!=null && out!=null && expected>0, "Bind frozen prediction metadata/context and exact population");Shared.setThreads(1);
		try{final CobalaminGenomicContext x=new CobalaminGenomicContext();x.loci(loci, expected);x.context(context);x.annotations();x.write(out);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Source binding      ----------------*/
	/*--------------------------------------------------------------*/

	private void loci(String path, int expected){
		read(path, NcrnaReviewUnion.LOCUS_HEADER, p->{
			final Locus l=new Locus(p);require(!byId.contains(l.id), "Duplicate review ID would merge distinct predictions");byId.put(l.id, l);loci.add(l);
			Group g=byContig.get(l.seqid);
			if(g==null){g=new Group(l.genome, l.seqid, l.length);byContig.put(l.seqid, g);}
			require(g.genome.equals(l.genome) && g.length==l.length, "One contig must retain one source owner and length");g.loci.add(l);l.group=g;
		});
		require(loci.size()==expected, "Frozen review population changed");
	}
	private void context(String path){
		read(path, NcrnaReviewJoin.CONTEXT_HEADER, p->{
			final Locus l=byId.get(p.parseString(0));
			require(l!=null && !l.bound && p.termEquals(l.genome, 1) && p.termEquals(l.assembly, 2) && p.termEquals(l.seqid, 3) && p.parseInt(4)==l.start
				&& p.parseInt(5)==l.stop && p.termEquals(Character.toString(l.strand), 6) && p.parseInt(7)==l.stop-l.start+1 && p.termEquals(l.fasta, 24),
				"Context must bind each exact original prediction and source assembly once");
			final String pathGff=p.parseString(22), hash=p.parseString(23);require(hash.matches("[0-9a-f]{20}"), "Source provenance must be sha80");
			l.bound=true;l.expectedFeatures=p.parseInt(11);l.expectedRegions=p.parseInt(12);
			final Group g=l.group;
			if(g.gff==null){g.gff=pathGff;g.hash=hash;}else{require(g.gff.equals(pathGff) && g.hash.equals(hash), "Conflicting GFF source bindings");}
			final String existing=files.get(pathGff);
			if(existing==null){files.put(pathGff, hash);paths.add(pathGff);}else{require(existing.equals(hash), "One GFF must have one bound sha80");}
		});
		for(Locus l:loci){require(l.bound, "Missing source context cannot become an intergenic result");}
	}
	private void annotations(){
		for(String path:paths){
			final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=0;String last=null;Group group=null;
			try{for(byte[] row; (row=in.nextLine())!=null; ){
				line++;if(row.length==0){continue;}if(row[0]=='#'){if(Tools.startsWith(row, "##FASTA")){break;}continue;}
				p.set(row);require(p.terms()==9, "Original GFF9 schema required");
				if(last==null || !p.termEquals(last, 0)){last=p.parseString(0);group=byContig.get(last);}
				if(group==null){continue;}require(group.gff.equals(path), "A source GFF may not borrow another assembly's contig");
				final int start=p.parseInt(3), stop=p.parseInt(4);final String strand=p.parseString(6);
				require(start>0 && start<=group.length && stop>=start && (long)stop-start+1<=group.length
					&& (strand.equals("+") || strand.equals("-") || strand.equals(".") || strand.equals("?")), "Invalid original feature geometry/strand");
				if(p.termEquals("region", 2)){
					final String circular=attribute(p.parseString(8), "Is_circular");
					if(circular!=null){
						require(start==1 && stop==group.length && (circular.equals("true") || circular.equals("false")), "Circular topology must describe the exact original contig");
						final boolean value=circular.equals("true");require(!group.topologyKnown || group.circular==value, "Conflicting source topology");
						group.topologyKnown=true;group.circular=value;
					}
					for(Locus l:group.loci){if(NcrnaReviewFeatureContext.overlap(l.start, l.stop, start, stop, group.length)>0){l.regions++;}}
					continue;
				}
				if(stop>group.length){group.wrapped=true;}
				final boolean cds=p.termEquals("CDS", 2);Feature f=null;
				for(Locus l:group.loci){
					final int overlap=NcrnaReviewFeatureContext.overlap(l.start, l.stop, start, stop, group.length);
					// Retain the bounded possible circular neighbor; choose only after the
					// complete GFF has supplied topology, irrespective of row ordering.
					final int gap=cds?downstreamGap(l.start, l.stop, l.strand, start, stop, group.length, true):-1;
					if(overlap==0 && (gap<0 || gap>300)){continue;}
					if(f==null){f=new Feature(p, line);}
					if(overlap>0){l.overlaps.add(new Overlap(f, overlap));}
					if(gap>=0 && gap<=300){l.candidates.add(f);}
				}
			}}finally{require(!in.close(), "Original GFF read failed: "+path);}
		}
		for(Locus l:loci){
			require(!l.group.wrapped || l.group.topologyKnown && l.group.circular, "Origin-spanning features require explicit circular topology");
			require(l.overlaps.size()==l.expectedFeatures && l.regions==l.expectedRegions,
				"New source pass must reproduce the accepted caller-overlap inventory: "+l.id);
			for(Feature f:l.candidates){
				final int gap=downstreamGap(l.start, l.stop, l.strand, f.start, f.stop, l.length, l.group.circular);
				if(gap>=0 && gap<=300){l.any.add(f, gap);if(f.strand==l.strand){l.same.add(f, gap);}}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String out){
		RrnaResourceIO.fresh(out, ".neighborhood.tsv", ".overlaps.tsv", ".neighbors.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".neighborhood.tsv"), overlaps=RrnaResourceIO.writer(out+".overlaps.tsv"), neighbors=RrnaResourceIO.writer(out+".neighbors.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			summary.println(HEADER);
			overlaps.println("review_id\tfeature_class\tstrand_relation\toverlap_bp\tfraction_of_call\tgff\tgff_sha80\tgff_line\ttype\tstart1\tstop1\tstrand\tgene\tproduct\toriginal_attributes");
			neighbors.println("review_id\tselection\tdistance_bp\ttied_nearest\tgff\tgff_sha80\tgff_line\tstart1\tstop1\tstrand\tgene\tproduct\tB12_name_hint\thint_rule\toriginal_attributes");
			for(Locus l:loci){
				int mask=0, parents=0, same=0, opposite=0;final int length=l.stop-l.start+1;
				for(Overlap o:l.overlaps){
					final Feature f=o.feature;final String relation=f.strand=='.' || f.strand=='?'?"UNKNOWN":f.strand==l.strand?"SAME":"OPPOSITE";
					if(f.parent){parents++;}else if(f.type.equals("CDS")){
						if(relation.equals("SAME")){mask|=1;same=Math.max(same, o.bases);}else if(relation.equals("OPPOSITE")){mask|=2;opposite=Math.max(opposite, o.bases);}else{mask|=4;}
					}else if(f.type.equals("tRNA")){mask|=8;}else if(f.type.equals("rRNA")){mask|=16;}else if(f.rna){mask|=32;}else{mask|=64;}
					overlaps.print(b.clear().append(l.id).tab().append(f.parent?"PARENT":f.type.equals("CDS")?"CDS":f.type.equals("tRNA")?"tRNA":f.type.equals("rRNA")?"rRNA":f.rna?"OTHER_NCRNA":"OTHER_FEATURE")
						.tab().append(relation).tab().append(o.bases).tab().append((double)o.bases/length, 8).tab().append(l.group.gff).tab().append(l.group.hash)
						.tab().append(f.line).tab().append(f.type).tab().append(f.start).tab().append(f.stop).tab().append(f.strand).tab().append(f.gene).tab().append(f.product).tab().append(f.attrs).nl());
				}
				b.clear().append(l.id).tab();classes(b, mask, parents);
				b.tab().append(same).tab().append((double)same/length, 8).tab().append(opposite).tab().append((double)opposite/length, 8).tab().append(same==length);
				neighborSummary(b, l.same, false);neighborSummary(b, l.any, true);summary.print(b.nl());
				neighborRows(neighbors, b, l, l.same, "SAME_STRAND");neighborRows(neighbors, b, l, l.any, "ANY_STRAND");
			}
		}finally{RrnaResourceIO.finish(summary, overlaps, neighbors);}
		System.err.println("COBALAMIN_GENOMIC_CONTEXT_PASS loci="+loci.size()+" gffs="+paths.size()+"; strict downstream0..300bp; names are hints only");
	}
	private static void neighborSummary(ByteBuilder b, Nearest n, boolean strand){
		b.tab();if(n.features.isEmpty()){b.append("NA\tfalse\tNA\tNA\t0");if(strand){b.append("\tNA");}return;}
		final Feature first=n.features.get(0);
		b.append(n.gap).tab().append(n.hint()).tab().append(first.gene).tab().append(first.product).tab().append(n.features.size());
		if(strand){b.tab().append(first.strand);}
	}
	private static void neighborRows(ByteStreamWriter out, ByteBuilder b, Locus l, Nearest n, String selection){
		for(Feature f:n.features){out.print(b.clear().append(l.id).tab().append(selection).tab().append(n.gap).tab().append(n.features.size())
			.tab().append(l.group.gff).tab().append(l.group.hash).tab().append(f.line).tab().append(f.start).tab().append(f.stop).tab().append(f.strand)
			.tab().append(f.gene).tab().append(f.product).tab().append(!f.rule.equals("NONE")).tab().append(f.rule).tab().append(f.attrs).nl());}
	}
	private static void classes(ByteBuilder b, int mask, int parents){
		if(mask==0){b.append(parents>0?"PARENT_ONLY":"INTERGENIC");return;}
		boolean first=true;for(int i=0; i<CLASSES.length; i++){if((mask&(1<<i))!=0){if(!first){b.append('+');}b.append(CLASSES[i]);first=false;}}
	}
	/** Strictly nonoverlapping CDS in transcript direction; adjacent gap is zero.
	 * A circular route is allowed only from explicit source topology. */
	public static int downstreamGap(int start, int stop, char strand, int geneStart, int geneStop, int length, boolean circular){
		require(start>0 && stop>=start && stop<=length && (strand=='+' || strand=='-') && geneStart>0 && geneStart<=length
			&& geneStop>=geneStart && (long)geneStop-geneStart+1<=length, "Downstream distances require valid original single-turn intervals");
		if(NcrnaReviewFeatureContext.overlap(start, stop, geneStart, geneStop, length)>0){return -1;}
		if(strand=='+'){
			if(geneStart>stop){return geneStart-stop-1;}return circular?length-stop+geneStart-1:-1;
		}
		final int end=geneStop>length?geneStop-length:geneStop;
		if(end<start){return start-end-1;}return circular?start+length-end-1:-1;
	}
	/** Literal requested gene-name families and B12-specific product terms only. */
	public static String hint(String gene, String product){
		final String g=gene.toLowerCase(Locale.ROOT), p=product.toLowerCase(Locale.ROOT);
		for(String prefix:GENES){if(g.startsWith(prefix)){
			final int n=prefix.length();
			if(prefix.equals("cob") || prefix.equals("cbi") || prefix.equals("cbl") || g.length()==n
				|| g.charAt(n)=='_' || g.charAt(n)=='-' || Character.isDigit(g.charAt(n))){return "gene:"+prefix;}
		}}
		if(p.contains("cobalamin")){return "product:cobalamin";}if(p.contains("cobinamide")){return "product:cobinamide";}
		if(p.contains("vitamin b12") || p.contains("vitamin%20b12") || p.contains("b12-dependent")){return "product:B12";}
		return "NONE";
	}
	private static void read(String path, String header, java.util.function.Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.UTF_8));final int terms=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);require(p.terms()==terms, "Frozen input schema changed");action.accept(p);}}
		finally{require(!in.close(), "Frozen input read failed");}
	}
	/** GFF values end at semicolons, not spaces; preserve raw product wording. */
	public static String attribute(String attrs, String key){
		require(attrs!=null && key!=null && !key.isEmpty(), "Read an exact GFF attribute key from intact source attributes");
		String value=null;final String prefix=key+'=';
		for(int from=0; from<attrs.length(); ){
			int end=attrs.indexOf(';', from);if(end<0){end=attrs.length();}
			if(end-from>=prefix.length() && attrs.startsWith(prefix, from)){
				require(value==null, "Repeated GFF attribute cannot silently choose a different product/gene: "+key);value=attrs.substring(from+prefix.length(), end);
			}
			from=end+1;
		}return value;
	}
	private static String field(String attrs, String key){final String value=attribute(attrs, key);return value==null?"NA":value;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Locus {
		Locus(LineParser1 p){
			id=p.parseString(0);genome=p.parseString(1);assembly=p.parseString(2);fasta=p.parseString(3);seqid=p.parseString(4);final String s=p.parseString(6);
			require(s.equals("+") || s.equals("-"), "Prediction strand must be explicit");strand=s.charAt(0);start=p.parseInt(7);stop=p.parseInt(8);length=p.parseInt(9);
			require(!id.isEmpty() && p.termEquals("cobalamin", 5) && p.termEquals("PREDICTION_UNRESOLVED", 13)
				&& start>0 && stop>=start && stop<=length, "Use valid frozen prediction intervals");
		}
		final String id, genome, assembly, fasta, seqid;final int start, stop, length;final char strand;Group group;boolean bound;int expectedFeatures, expectedRegions, regions;
		final ArrayList<Overlap> overlaps=new ArrayList<Overlap>();final ArrayList<Feature> candidates=new ArrayList<Feature>();final Nearest same=new Nearest(), any=new Nearest();
	}
	private static final class Group {
		Group(String g, String i, int n){genome=g;id=i;length=n;}final String genome, id;final int length;String gff, hash;boolean topologyKnown, circular, wrapped;
		final ArrayList<Locus> loci=new ArrayList<Locus>();
	}
	private static final class Feature {
		Feature(LineParser1 p, long n){
			line=n;type=p.parseString(2);start=p.parseInt(3);stop=p.parseInt(4);strand=p.parseString(6).charAt(0);attrs=p.parseString(8);
			gene=field(attrs, "gene");product=field(attrs, "product");rule=hint(gene, product);
			parent=type.equals("gene") || type.equals("pseudogene") || type.equals("mRNA") || type.equals("transcript");
			rna=type.endsWith("RNA") || type.endsWith("_RNA") || type.equals("riboswitch") || attribute(attrs, "ncRNA_class")!=null;
		}
		final long line;final String type, attrs, gene, product, rule;final char strand;final int start, stop;final boolean parent, rna;
	}
	private static final class Overlap {Overlap(Feature f, int n){feature=f;bases=n;}final Feature feature;final int bases;}
	private static final class Nearest {
		void add(Feature f, int distance){
			assert(distance>=0 && distance<=300):"Only the requested downstream interval may enter nearest-CDS selection";
			if(distance>gap){return;}if(distance<gap){features.clear();gap=distance;}features.add(f);
		}
		boolean hint(){for(Feature f:features){if(!f.rule.equals("NONE")){return true;}}return false;}
		int gap=Integer.MAX_VALUE;final ArrayList<Feature> features=new ArrayList<Feature>();
	}
	private final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ObjectMap<String,Group> byContig=new ObjectMap<String,Group>(String.class, Group.class);
	private final ObjectMap<String,String> files=new ObjectMap<String,String>(String.class, String.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<String> paths=new ArrayList<String>();
	private static final String[] GENES={"btub", "btuf", "cob", "cbi", "mete", "meth", "nrdj", "pduq", "euta", "cbl", "succ"};
	private static final String[] CLASSES={"CDS_SAME", "CDS_OPPOSITE", "CDS_UNKNOWN", "tRNA", "rRNA", "OTHER_NCRNA", "OTHER_FEATURE"};
	public static final String HEADER="review_id\toverlap_class\tmax_CDS_same_overlap_bp\tmax_CDS_same_fraction\tmax_CDS_opposite_overlap_bp\tmax_CDS_opposite_fraction\tCDS_same_contains_call"
		+"\tnearest_downstream_same_distance_bp\tnearest_downstream_same_B12_hint_any_tied\tnearest_downstream_same_gene_first\tnearest_downstream_same_product_first\tnearest_downstream_same_ties"
		+"\tnearest_downstream_any_distance_bp\tnearest_downstream_any_B12_hint_any_tied\tnearest_downstream_any_gene_first\tnearest_downstream_any_product_first\tnearest_downstream_any_ties\tnearest_downstream_any_strand_first";
}
