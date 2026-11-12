package prok;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AccessDeniedException;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;

import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaTextIO.Input;
import structures.ByteBuilder;
import structures.IntList;

/** Retains source GFF evidence and RNA-related structure without admitting truth.
 * One source file owns its ID/Parent namespace. Original rows and line numbers
 * survive; multipart IDs never become a synthesized interval. Unknown feature
 * types remain reviewable, not presumed non-coding. Generic coding/region rows
 * survive only when directly tagged or needed by a selected structural branch.
 * Working memory covers one GFF. Unselected raw rows and the temporary graph are
 * released before return; the immutable source remains the exclusion audit trail.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaHistoricalEvidence{

	/*--------------------------------------------------------------*/
	/*----------------          Source input        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads a declared source; unavailable/unresolved paths are explicit evidence states.
	 * Malformed GFF structure or attributes fail loudly rather than producing a partial source.
	 */
	public static Source read(final String file, final Observation observed){
		require(file!=null && observed!=null, "Historical evidence requires a source claim and input-scoped contigs");
		final Source source=new Source(file, observed);
		if(file.equals("NA")){
			source.state="MISSING_GFF";
			return source;
		}
		final Path path=Paths.get(file);
		if(!path.isAbsolute()){
			source.state="UNRESOLVED_GFF_PATH";
			return source;
		}
		try{
			final BasicFileAttributes info=Files.readAttributes(path, BasicFileAttributes.class);
			if(!info.isRegularFile()){
				source.state="NOT_REGULAR_GFF";
				return source;
			}
			if(!Files.isReadable(path)){
				source.state="UNREADABLE_GFF";
				return source;
			}
		}catch(NoSuchFileException missing){
			source.state="MISSING_GFF";
			return source;
		}catch(AccessDeniedException denied){
			source.state="UNREADABLE_GFF";
			return source;
		}catch(IOException unresolved){
			source.state="UNRESOLVED_GFF_IO";
			return source;
		}
		final Graph graph=new Graph(source);
		final LineParser1 p=new LineParser1('\t');
		final Attributes attributes=new Attributes();
		boolean inline=false;
		String encodedTarget=null, target=null;
		try(Input input=new Input(file)){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				source.lines++;
				if(inline){
					source.inlineFasta++;
					continue;
				}
				if(row.length==0){
					source.blank++;
					continue;
				}
				if(row[0]=='#'){
					source.comments++;
					inline=EukOrganellePanelInventory.starts(row, "##FASTA") && (row.length==7 || row[7]<=32);
					continue;
				}
				p.set(row);
				if(p.terms()!=9){throw new IllegalArgumentException("Expected nine GFF fields at "+file+":"+source.lines);}
				if(encodedTarget==null || !p.termEquals(encodedTarget, 0)){
					encodedTarget=text(p, 0);
					p.setBounds(0);
					target=attributes.decode(row, p.a(), p.b());
				}
				attributes.read(p);
				final String type=type(p);
				int reason=isRna(type) || attributes.rnaHint ? DIRECT_RNA : 0;
				if(attributes.rfam.length>0){reason|=DIRECT_RFAM;}
				if(attributes.unknownBiotype || (!isRna(type) && exclusion(type)==null)){reason|=DIRECT_UNKNOWN;}
				final long start=EukPanelGffValidation.coordinate(p, 3), stop=EukPanelGffValidation.coordinate(p, 4);
				p.setBounds(6);
				final char strand=p.b()==p.a()+1 ? (char)row[p.a()] : '?';
				final Feature feature=new Feature(source.lines, row, target, type, start, stop, strand,
					attributes.id, attributes.parents, attributes.rfam, reason, observed.length(target));
				graph.add(feature);
			}
		}catch(IllegalArgumentException failure){
			throw new IllegalArgumentException(file+":"+source.lines+": "+failure.getMessage(), failure);
		}
		graph.retain();
		source.state=graph.rows.isEmpty() ? "EMPTY_GFF" : "PARSED_GFF";
		require(source.lines==source.comments+source.blank+source.inlineFasta+source.features,
			"Every source line must remain accounted for: "+file);
		require(source.features==source.retained.size()+source.excluded,
			"Retained rows and explicit exclusions must conserve source features: "+file);
		return source;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Structural graph       ----------------*/
	/*--------------------------------------------------------------*/

	/** Temporary ID graph; descendant and ancestor traversal have separate selection state. */
	private static final class Graph{
		Graph(final Source source_){source=source_;}

		void add(final Feature feature){
			assert(rows.size()==source.features) : "Each parsed GFF row enters the graph once";
			final int row=rows.size();
			final Node node=feature.id==null ? node(null) : identified(feature.id);
			if(node.rows.size>0){
				final Feature first=rows.get(node.rows.get(0));
				node.incompatible|=!first.target.equals(feature.target) || !first.type.equals(feature.type) || first.strand!=feature.strand;
			}
			rows.add(feature);
			rowNodes.add(node.index);
			node.rows.add(row);
			node.direct|=feature.selection!=0;
			for(final String parent:feature.parents){
				final Node ancestor=identified(parent);
				if(!node.parents.contains(ancestor.index)){
					node.parents.add(ancestor.index);
					node.parentRows.add(row);
					ancestor.children.add(node.index);
				}
			}
			source.features++;
		}

		private Node identified(final String id){
			Node node=ids.get(id);
			if(node==null){
				node=node(id);
				ids.put(id, node);
			}
			return node;
		}

		private Node node(final String id){
			final Node node=new Node(nodes.size(), id);
			nodes.add(node);
			return node;
		}

		void retain(){
			assert(rows.size()==source.features) : "Selection operates on the complete source graph";
			final IntList queue=new IntList();
			for(final Node node:nodes){
				if(node.direct){
					node.descendant=true;
					queue.add(node.index);
				}
			}
			for(int i=0; i<queue.size; i++){
				final Node current=nodes.get(queue.get(i));
				for(int j=0; j<current.children.size; j++){
					final Node child=nodes.get(current.children.get(j));
					if(!child.descendant){
						child.descendant=true;
						queue.add(child.index);
					}
				}
			}
			// Only the selected branch seeds this upward walk. Context ancestors do
			// not seed another downward walk into sibling CDS or transcript branches.
			for(int i=0; i<queue.size; i++){
				final Node current=nodes.get(queue.get(i));
				for(int j=0; j<current.parents.size; j++){
					final Node parent=nodes.get(current.parents.get(j));
					if(!parent.ancestor){
						parent.ancestor=true;
						if(!parent.descendant){queue.add(parent.index);}
					}
				}
			}
			for(final Node node:nodes){
				if(node.rows.size==0 && (node.descendant || node.ancestor)){source.missingParents.add(node.id);}
			}
			cycles();
			for(int i=0; i<rows.size(); i++){
				final Feature feature=rows.get(i);
				final Node node=nodes.get(rowNodes.get(i));
				if(node.descendant || node.ancestor){
					if(node.descendant && feature.selection==0){feature.selection|=DESCENDANT;}
					if(node.ancestor){feature.selection|=ANCESTOR;}
					feature.incompatibleId=node.incompatible;
					source.retained.add(feature);
				}else{
					final String reason=exclusion(feature.type);
					if(reason==null){throw new IllegalStateException("An unknown source type must remain reviewable: "+feature.type);}
					final int count=source.exclusions.contains(reason) ? source.exclusions.get(reason) : 0;
					source.exclusions.put(reason, count+1);
					source.excluded++;
				}
			}
		}

		/** Records directed DFS back-edge witnesses, not a claim to enumerate every cycle or SCC. */
		private void cycles(){
			final IntList stack=new IntList(), next=new IntList();
			for(final Node root:nodes){
				if(root.color!=0){continue;}
				stack.add(root.index);
				next.add(0);
				root.color=1;
				while(stack.size>0){
					assert(stack.size==next.size) : "Each DFS frame owns one next-parent position";
					final int depth=stack.size-1;
					final Node node=nodes.get(stack.get(depth));
					final int edge=next.get(depth);
					if(edge==node.parents.size){
						node.color=2;
						stack.pop();
						next.pop();
						continue;
					}
					next.set(depth, edge+1);
					final Node parent=nodes.get(node.parents.get(edge));
					if(parent.color==0){
						parent.color=1;
						stack.add(parent.index);
						next.add(0);
					}else if(parent.color==1){
						source.cycleWitnesses.add(new Relation(node.id, parent.id, rows.get(node.parentRows.get(edge)).line));
					}
				}
			}
		}

		final Source source;
		final ArrayList<Feature> rows=new ArrayList<Feature>();
		final IntList rowNodes=new IntList();
		final ArrayList<Node> nodes=new ArrayList<Node>();
		final ObjectMap<String,Node> ids=new ObjectMap<String,Node>(String.class, Node.class);
	}

	private static final class Node{
		Node(final int index_, final String id_){
			index=index_;
			id=id_;
		}
		final int index;
		final String id;
		// Most IDs occupy one row and one Parent edge; the default capacity is 256.
		final IntList rows=new IntList(1), parents=new IntList(1), parentRows=new IntList(1), children=new IntList(1);
		boolean direct, descendant, ancestor, incompatible;
		int color;
	}

	/*--------------------------------------------------------------*/
	/*----------------        GFF attributes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reused byte scanner; only retained identifiers/relationships are materialized. */
	static final class Attributes{
		void read(final LineParser1 p){
			assert(p.terms()==9) : "Attribute parsing requires the validated ninth GFF field";
			id=null;
			parents=EMPTY;
			rfam=EMPTY;
			rnaHint=false;
			unknownBiotype=false;
			boolean parentSeen=false, idSeen=false;
			ArrayList<String> models=null;
			p.setBounds(8);
			final byte[] row=p.line();
			final int stop=p.b();
			if(stop==p.a()+1 && row[p.a()]=='.'){return;}
			for(int from=p.a(); from<stop; ){
				int end=from;
				while(end<stop && row[end]!=';'){end++;}
				int equal=from;
				while(equal<end && row[equal]!='='){equal++;}
				require(end==from || (equal>from && equal<end), "Malformed GFF attribute: expected name=value");
				if(equal<end){
					if(same(row, from, equal, "ID")){
						require(!idSeen, "Duplicate GFF ID attribute");
						id=decode(row, equal+1, end);
						require(!id.isEmpty(), "Empty GFF ID attribute");
						idSeen=true;
					}else if(same(row, from, equal, "Parent")){
						require(!parentSeen, "Duplicate GFF Parent attribute");
						parents=values(row, equal+1, end);
						parentSeen=true;
					}else if(same(row, from, equal, "Dbxref")){
						for(int begin=equal+1; begin<=end; ){
							int comma=begin;
							while(comma<end && row[comma]!=','){comma++;}
							if(comma-begin>=5 && same(row, begin, begin+4, "RFAM")){
								final String value=decode(row, begin, comma);
								if(value.startsWith("RFAM:")){
									require(value.length()>5, "Empty RFAM cross-reference");
									if(models==null){models=new ArrayList<String>(1);}
									models.add(value.substring(5));
								}
							}
							begin=comma+1;
						}
					}else if(same(row, from, equal, "ncRNA_class")){
						rnaHint=true;
					}else if(same(row, from, equal, "gene_biotype") || same(row, from, equal, "transcript_biotype")){
						final int hint=hint(row, equal+1, end);
						rnaHint|=hint==DIRECT_RNA;
						unknownBiotype|=hint==DIRECT_UNKNOWN;
					}else if(same(row, from, equal, "gbkey")){
						// gbkey=Gene is a structural annotation, not an unknown biotype.
						rnaHint|=hint(row, equal+1, end)==DIRECT_RNA;
					}
				}
				from=end+1;
			}
			if(models!=null){rfam=models.toArray(new String[models.size()]);}
		}

		/** Common ASCII claims stay on bytes; escaped/Unicode claims use the strict decoder.
		 * Returns a direct RNA/unknown flag, or zero for explicit protein-coding/mRNA.
		 */
		private int hint(final byte[] row, final int start, final int stop){
			assert(start>=0 && start<=stop && stop<=row.length) : "Hint classification stays inside its parsed attribute value";
			for(int i=start; i<stop; i++){
				if(row[i]=='%' || row[i]<0){
					final String value=decode(row, start, stop);
					if(isRna(value)){return DIRECT_RNA;}
					return value.equals("protein_coding") || value.equals("mRNA") ? 0 : DIRECT_UNKNOWN;
				}
			}
			if(same(row, start, stop, "protein_coding") || same(row, start, stop, "mRNA")){return 0;}
			return stop-start>=3 && same(row, stop-3, stop, "RNA") ? DIRECT_RNA : DIRECT_UNKNOWN;
		}

		private String[] values(final byte[] row, final int start, final int stop){
			int count=1;
			for(int i=start; i<stop; i++){if(row[i]==','){count++;}}
			final String[] result=new String[count];
			int index=0;
			for(int from=start; from<=stop; ){
				int end=from;
				while(end<stop && row[end]!=','){end++;}
				require(end>from, "Empty GFF Parent identifier");
				result[index++]=decode(row, from, end);
				from=end+1;
			}
			assert(index==count) : "Comma splitting precedes percent decoding, so encoded commas stay inside identifiers";
			return result;
		}

		String decode(final byte[] row, final int start, final int stop){
			assert(start>=0 && start<=stop && stop<=row.length) : "Decoding is bounded by one validated field or attribute value";
			buffer.clear();
			for(int i=start; i<stop; i++){
				if(row[i]=='%'){
					require(i+2<stop, "Incomplete GFF percent escape");
					final int a=hex(row[++i]), b=hex(row[++i]);
					require(a>=0 && b>=0, "Invalid GFF percent escape");
					buffer.append((byte)((a<<4)|b));
				}else{buffer.append(row[i]);}
			}
			final String decoded=new String(buffer.array, 0, buffer.length(), StandardCharsets.UTF_8);
			require(decoded.indexOf('\uFFFD')<0, "GFF identifier loses identity through invalid UTF-8 or a replacement character");
			return decoded;
		}

		String id;
		String[] parents=EMPTY, rfam=EMPTY;
		boolean rnaHint, unknownBiotype;
		final ByteBuilder buffer=new ByteBuilder();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Results            ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Source{
		private Source(final String file_, final Observation observed_){
			file=file_;
			observed=observed_;
		}
		public int size(){return retained.size();}
		public Feature get(final int index){return retained.get(index);}
		public String state(){return state;}
		public long featureRows(){return features;}
		public long excludedRows(){return excluded;}
		public int excluded(final String reason){return exclusions.contains(reason) ? exclusions.get(reason) : 0;}
		public int missingParents(){return missingParents.size();}
		public String missingParent(final int index){return missingParents.get(index);}
		public int cycles(){return cycleWitnesses.size();}
		public Relation cycle(final int index){return cycleWitnesses.get(index);}
		public final String file;
		public final Observation observed;
		private String state;
		private long lines, comments, blank, inlineFasta, features, excluded;
		private final ArrayList<Feature> retained=new ArrayList<Feature>();
		private final ArrayList<String> missingParents=new ArrayList<String>();
		private final ArrayList<Relation> cycleWitnesses=new ArrayList<Relation>();
		private final ObjectIntMap<String> exclusions=new ObjectIntMap<String>(String.class);
	}

	public static final class Feature{
		private Feature(final long line_, final byte[] raw_, final String target_, final String type_, final long start_,
			final long stop_, final char strand_, final String id_, final String[] parents_, final String[] rfam_, final int reason, final int length){
			line=line_;
			raw=raw_;
			target=target_;
			type=type_;
			start=start_;
			stop=stop_;
			strand=strand_;
			id=id_;
			parents=parents_;
			rfam=rfam_;
			selection=reason;
			range=start<1 || stop<start ? "INVALID_COORDINATES" : length<0 ? "UNRESOLVED_LITERAL_TARGET"
				: stop>length ? "OUT_OF_RANGE" : "VALID_INTERVAL";
		}
		public int selection(){return selection;}
		public boolean incompatibleId(){return incompatibleId;}
		public final long line, start, stop;
		/** Original bytes; callers must not mutate returned records or identifier arrays. */
		public final byte[] raw;
		public final String target, type, id, range;
		public final char strand;
		public final String[] parents, rfam;
		private int selection;
		private boolean incompatibleId;
	}

	/** One actual directed Parent back-edge witness, retaining its original source line. */
	public static final class Relation{
		private Relation(final String child_, final String parent_, final long line_){
			child=child_;
			parent=parent_;
			line=line_;
		}
		public final String child, parent;
		public final long line;
	}

	static String type(final LineParser1 p){
		for(final String type:COMMON_TYPES){if(p.termEquals(type, 2)){return type;}}
		return text(p, 2);
	}
	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static boolean isRna(final String type){return type.endsWith("RNA") && !type.equals("mRNA");}
	private static String exclusion(final String type){
		if(type.equals("CDS") || type.equals("mRNA")){return "UNSELECTED_CODING";}
		if(type.equals("gene") || type.equals("exon")){return "UNSELECTED_GENERIC_STRUCTURE";}
		return type.equals("region") ? "UNSELECTED_REGION" : null;
	}
	private static boolean same(final byte[] row, final int start, final int stop, final String text){
		if(stop-start!=text.length()){return false;}
		for(int i=0; i<text.length(); i++){if(row[start+i]!=text.charAt(i)){return false;}}
		return true;
	}
	private static int hex(final byte value){
		if(value>='0' && value<='9'){return value-'0';}
		if(value>='A' && value<='F'){return value-'A'+10;}
		return value>='a' && value<='f' ? value-'a'+10 : -1;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private NcrnaHistoricalEvidence(){}
	public static final int DIRECT_RNA=1, DIRECT_RFAM=2, DIRECT_UNKNOWN=4, DESCENDANT=8, ANCESTOR=16;
	private static final String[] EMPTY=new String[0];
	private static final String[] COMMON_TYPES={"CDS", "gene", "mRNA", "exon", "region", "tRNA", "rRNA", "ncRNA", "tmRNA"};
}
