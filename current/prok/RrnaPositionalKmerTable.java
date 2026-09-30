package prok;

import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import parse.Parse;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Per-consensus, per-end log enrichment of the anchor k-mer against the other
 * 24 positions. Offsets are k-mer START coordinates relative to the true first
 * or last gene base in sense-oriented flanked reads, not alignment columns.
 * Uses native streaming/count arrays; no network training or caller changes.
 * @author Ganyu
 */
public final class RrnaPositionalKmerTable {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaPositionalKmerTable.class,false);
		String mode=null,in=null,model=null,prefix=null,table=null,out=null;int k=-1,r=-1,start=Integer.MIN_VALUE,stop=Integer.MIN_VALUE;long expected=-1;double alpha=Double.NaN;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected key=value");final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(key.equals("mode")){mode=v;}else if(key.equals("in")){in=v;}else if(key.equals("model")){model=v;}else if(key.equals("outprefix")){prefix=v;}else if(key.equals("table")){table=v;}else if(key.equals("out")){out=v;}
			else if(key.equals("k")){k=Integer.parseInt(v);}else if(key.equals("r")||key.equals("radius")){r=Integer.parseInt(v);}else if(key.equals("startoffset")){start=Integer.parseInt(v);}else if(key.equals("stopoffset")){stop=Integer.parseInt(v);}
			else if(key.equals("expected")){expected=Long.parseLong(v);}else if(key.equals("pseudocount")){alpha=Double.parseDouble(v);}else{throw new IllegalArgumentException("Unknown option: "+key);}}
		require(in!=null&&model!=null&&expected>0,"Declare input, assigned consensus and exact record count");
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;
		if("build".equals(mode)){require(prefix!=null&&table==null&&out==null&&r==RADIUS&&start!=Integer.MIN_VALUE&&stop!=Integer.MIN_VALUE,"Build requires outprefix, r=12 and both signed k-mer-start offsets; table/out belong to score mode");build(in,model,prefix,expected,k,start,stop,alpha);}
		else if("score".equals(mode)){require(table!=null&&out!=null&&prefix==null&&k==-1&&r==-1&&start==Integer.MIN_VALUE&&stop==Integer.MIN_VALUE&&Double.isNaN(alpha),"Score takes model/k/geometry/smoothing from its saved table; do not supply ignored build overrides");score(in,load(table,model),out,expected);}
		else{throw new IllegalArgumentException("mode must be build or score");}
	}
	/** Input is one explicitly assigned consensus partition. An optional model=
	 * header must agree; callers must never pool differently assigned partitions. */
	static void build(String in,String model,String prefix,long expected,int k,int start,int stop,double alpha){
		final Table five=new Table(model,"5prime",k,start,alpha),three=new Table(model,"3prime",k,stop,alpha);
		RrnaResourceIO.fresh(prefix,".5prime.tsv.gz",".3prime.tsv.gz");final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(in,read->{checkModel(read.id,model);require(seen.add(RrnaResourceIO.first(read.id)),"Duplicate training record would alter table weights");
			final int left=Integer.parseInt(RrnaResourceIO.field(read.id,"lflank")),right=Integer.parseInt(RrnaResourceIO.field(read.id,"rflank"));
			require(read.bases!=null&&left>=0&&right>=0&&(long)left+right<read.length(),"Sense-oriented flanks must leave a nonempty true gene");
			five.add(read.bases,left);three.add(read.bases,read.length()-right-1);
		});require(five.records==expected&&three.records==expected,"Assigned training partition count changed");five.finish();three.finish();five.write(prefix+".5prime.tsv.gz");three.write(prefix+".3prime.tsv.gz");
		System.err.println("RRNA_POSITIONAL_TABLE_BUILD_PASS model="+model+" records="+expected+" k="+k+" radius=12 pseudocount="+alpha);
	}
	static void checkModel(String header,String model){final String assignment=RrnaResourceIO.optionalField(header,"model");require(assignment==null||assignment.equals(model),"Header consensus assignment differs from declared input partition");}
	public static final class Table{
		public Table(String model_,String end_,int k_,int offset_,double alpha_){
			require(model_!=null&&model_.matches("[A-Za-z0-9_.:-]+")&&(end_.equals("5prime")||end_.equals("3prime")),"One explicit consensus and gene end per table");
			require(k_>=6&&k_<=9&&Double.isFinite(alpha_)&&alpha_>0,"Recipe anchor search uses k6..9 and finite positive pseudocount");
			model=model_;end=end_;k=k_;offset=offset_;alpha=alpha_;space=1<<(2*k);counts=new long[SITES][space];
		}
		public void add(byte[] bases,int trueEnd){
			require(!finished&&bases!=null&&trueEnd>=0&&trueEnd<bases.length,"Training requires a real inclusive endpoint before finalization");records++;
			for(int site=0;site<SITES;site++){final long pos=(long)trueEnd+offset+site-RADIUS;
				if(pos<0||pos+k>bases.length){clipped[site]++;continue;}final int word=encode(bases,(int)pos,k);
				if(word<0){ambiguous[site]++;continue;}counts[site][word]++;valid[site]++;
			}
		}
		public void finish(){
			require(!finished&&records>0,"Finalize a nonempty training population once");
			for(int site=0;site<SITES;site++){long sum=0;for(long n:counts[site]){require(n>=0,"Count overflow must not produce a valid table");sum=Math.addExact(sum,n);}
				require(sum==valid[site]&&valid[site]+clipped[site]+ambiguous[site]==records,"Every training record has one counted, clipped or ambiguous observation per offset");
				require(valid[site]>0,"Every scored offset needs observed training bases, not only a pseudocount prior");
			}scores=new double[space];for(int word=0;word<space;word++){scores[word]=computeScore(word);}finished=true;
		}
		/** Add alpha to EVERY possible key at EVERY offset, normalize each offset
		 * separately, then average the 24 non-anchor frequencies (anchor excluded).
		 * Separate denominators matter when clipped or ambiguous sites differ. */
		public double scoreKey(int word){
			require(finished&&word>=0&&word<space,"Scoring requires a finalized table and valid encoded k-mer");
			return scores[word];
		}
		private double computeScore(int word){
			assert(word>=0&&word<space):"Finalize exactly one lookup score per encoded key";
			final double prior=alpha*space,anchor=(counts[RADIUS][word]+alpha)/(valid[RADIUS]+prior);double background=0;
			for(int site=0;site<SITES;site++){if(site!=RADIUS){background+=(counts[site][word]+alpha)/(valid[site]+prior);}}
			final double value=Math.log(anchor/(background/(SITES-1)));require(Double.isFinite(value),"Finite positive smoothing must produce a finite log ratio");return value;
		}
		/** anchorStart is the ZERO-BASED k-mer start at the center of the candidate
		 * window, already shifted by this table's endpoint-relative anchor offset.
		 * Invalid/clipped k-mers are NaN, never a fabricated neutral score. */
		public void scoreWindow(byte[] bases,int anchorStart,float[] dest){
			require(finished&&bases!=null&&dest!=null&&dest.length==SITES,"Emit exactly25 ordered candidate-position scores from a finalized table");
			for(int site=0;site<SITES;site++){final long pos=(long)anchorStart+site-RADIUS;final int word=pos<0||pos+k>bases.length?-1:encode(bases,(int)pos,k);dest[site]=word<0?Float.NaN:(float)scoreKey(word);}
		}
		public void write(String path){
			require(finished,"Only validated completed tables are written");RrnaResourceIO.fresh(path,"");final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
			try{out.println(SCHEMA);out.println(META_HEADER);out.print(b.clear().append(model).tab().append(end).tab().append(k).tab().append(RADIUS).tab().append(offset).tab().append(Double.toString(alpha)).tab().append(records).nl());
				writeTotals(out,b,"valid",valid);writeTotals(out,b,"clipped",clipped);writeTotals(out,b,"ambiguous",ambiguous);out.println(rowHeader());
				for(int word=0;word<space;word++){long other=0;for(int site=0;site<SITES;site++){if(site!=RADIUS){other+=counts[site][word];}}if(other+counts[RADIUS][word]==0){continue;}
					b.clear();appendKmer(b,word,k);b.tab().append(counts[RADIUS][word]).tab().append(other).tab().append(Double.toString(scoreKey(word)));
					for(int site=0;site<SITES;site++){b.tab().append(counts[site][word]);}out.print(b.nl());
				}
			}finally{RrnaResourceIO.finish(out);}
		}
		public final String model,end;public final int k,offset,space;public final double alpha;public long records;
		final long[][] counts;final long[] valid=new long[SITES],clipped=new long[SITES],ambiguous=new long[SITES];double[] scores;boolean finished;
	}
	public static Table load(String path,String expectedModel){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,SCHEMA,path);RrnaResourceIO.header(in,META_HEADER,path);p.set(in.nextLine());require(p.terms()==7&&p.termEquals(expectedModel,0)&&p.parseInt(3)==RADIUS,"Saved model and radius must match inference");
			final Table table=new Table(p.parseString(0),p.parseString(1),p.parseInt(2),p.parseInt(4),savedDouble(p,5));table.records=p.parseLong(6);
			readTotals(in,p,"valid",table.valid);readTotals(in,p,"clipped",table.clipped);readTotals(in,p,"ambiguous",table.ambiguous);RrnaResourceIO.header(in,rowHeader(),path);
			final boolean[] seen=new boolean[table.space];final double[] saved=new double[table.space];
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==SITES+4&&p.termStartsWith("",0)&&p.b()-p.a()==table.k,"Table rows must retain k-mer and all25 positional counts");final int word=encode(row,0,table.k);require(word>=0&&!seen[word],"Invalid/duplicate table key");seen[word]=true;long other=0;
				final long anchor=p.parseLong(1),background=p.parseLong(2);saved[word]=savedDouble(p,3);
				for(int site=0;site<SITES;site++){final long n=p.parseLong(site+4);require(n>=0,"Negative table count");table.counts[site][word]=n;if(site!=RADIUS){other=Math.addExact(other,n);}}
				require(anchor==table.counts[RADIUS][word]&&background==other&&anchor+other>0,"Anchor/background columns must match the saved positional counts");
			}table.finish();
			for(int word=0;word<table.space;word++){if(seen[word]){require(Double.isFinite(saved[word])&&Math.abs(saved[word]-table.scoreKey(word))<1e-10,"Saved log score disagrees with its counts/denominators/pseudocount");}}
			return table;
		}finally{require(!in.close(),"Score table read failed");}
	}
	static void score(String input,Table table,String path,long expected){
		RrnaResourceIO.fresh(path,"");final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();final float[] values=new float[SITES];final long[] rows={0};final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{b.append("record_id\tmodel\tend\tanchor_start\tvalid_positions");for(int site=0;site<SITES;site++){b.tab().append("score_").append(site-RADIUS);}out.print(b.nl());
			RrnaResourceIO.stream(input,read->{checkModel(read.id,table.model);final String id=RrnaResourceIO.first(read.id);require(seen.add(id),"Candidate IDs must be unique");final int anchor=Integer.parseInt(RrnaResourceIO.field(read.id,"anchor"));table.scoreWindow(read.bases,anchor,values);int valid=0;for(float v:values){if(Float.isFinite(v)){valid++;}}
				b.clear().append(id).tab().append(table.model).tab().append(table.end).tab().append(anchor).tab().append(valid);for(float value:values){b.tab();if(Float.isNaN(value)){b.append("NA");}else{b.appendSlow(value);}}out.print(b.nl());rows[0]++;
			});require(rows[0]==expected,"Candidate count changed; every row requires a25-position output");
		}finally{RrnaResourceIO.finish(out);}System.err.println("RRNA_POSITIONAL_TABLE_SCORE_PASS model="+table.model+" end="+table.end+" rows="+rows[0]);
	}
	static void writeTotals(ByteStreamWriter out,ByteBuilder b,String name,long[] counts){assert(counts.length==SITES):"One denominator per candidate offset";b.clear().append('#').append(name);for(long n:counts){b.tab().append(n);}out.print(b.nl());}
	/** Double.toString in write() can emit exponents. Keep native parsing for
	 * ordinary numbers and use the native slow fallback only at exponent fields. */
	static double savedDouble(LineParser1 p,int term){
		p.setBounds(term);final byte[] row=p.line();
		assert(p.a()<p.b()) : "Saved table numeric fields must be nonempty";
		for(int i=p.a();i<p.b();i++){if(row[i]=='E'||row[i]=='e'){return Parse.parseDoubleSlow(row,p.a(),p.b());}}
		return p.parseDouble(term);
	}
	static void readTotals(ByteFile in,LineParser1 p,String name,long[] counts){p.set(in.nextLine());require(p.terms()==SITES+1&&p.termEquals("#"+name,0),"Missing per-offset observation denominators: "+name);for(int i=0;i<SITES;i++){counts[i]=p.parseLong(i+1);require(counts[i]>=0,"Negative observation denominator");}}
	static String rowHeader(){final ByteBuilder b=new ByteBuilder().append("kmer\tanchor_count\tother_count\tlog_score");for(int site=0;site<SITES;site++){b.tab().append("count_").append(site-RADIUS);}return b.toString();}
	static int encode(byte[] bases,int from,int k){assert(from>=0&&(long)from+k<=bases.length):"Bounds are checked before encoding a candidate position";int word=0;for(int i=from;i<from+k;i++){final int n=RrnaResourceIO.base(bases[i]);if(n<0){return -1;}word=(word<<2)|n;}return word;}
	static void appendKmer(ByteBuilder b,int word,int k){for(int shift=2*(k-1);shift>=0;shift-=2){b.append(BASES[(word>>>shift)&3]);}}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final int RADIUS=12,SITES=25;
	static final String SCHEMA="#RRNA_ANCHOR_LOG_RATIO_V1",META_HEADER="model\tend\tk\tradius\tanchor_start_offset\tpseudocount_per_key_per_offset\trecords";
	static final byte[] BASES={'A','C','G','T'};
}
