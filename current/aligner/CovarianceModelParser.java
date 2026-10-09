package aligner;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import static aligner.CovarianceModel.*;

/** Strict single-model RNA INFERNAL1/a reader; filter HMM is explicitly out of scope.
 * Format/normalization: Infernal1.1.5 cm_file.c read_asc_1p1_cm/ascii2prob,
 * cm.c CMRenormalize/CMLogoddsify; consensus numbering: display.c CreateEmitMap.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelParser {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CovarianceModelParser.class, false);
		final Parser parser=new Parser();parser.out1="stdout.txt";
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			require(parser.parse(arg, key, value), "Unknown argument: "+key);}
		require(parser.in1!=null && parser.out1!=null && !parser.in1.equals(parser.out1), "Distinct in=CM and out=summary required");
		Shared.setThreads(1);final CovarianceModel m=read(parser.in1);
		final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(parser.out1, FileFormat.TEXT, null, true, parser.overwrite, false, false));
		w.start();final ByteBuilder b=new ByteBuilder();
		w.println("name\taccession\tstates\tnodes\tclen\twindow\tbps\tbifs\tD\tMP\tML\tMR\tIL\tIR\tS\tE\tB\tEL");
		b.append(m.name).tab().append(m.accession==null ? "NA" : m.accession).tab().append(m.states()).tab().append(m.nodes()).tab().append(m.clen).tab().append(m.window);
		b.tab().append(m.stateCount("MP")).tab().append(m.stateCount("B"));
		for(String state:STATE_NAMES){b.tab().append(m.stateCount(state));}w.print(b.nl());
		require(!w.poisonAndWait(), "CM summary output failed");Shared.closeStream(pp.outstream);
	}

	public static CovarianceModel read(String path){
		assert(path!=null):"CM reading requires a named input for diagnostics";
		final ByteFile in=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		try{return new CovarianceModelParser(path, in).parse();}
		finally{require(!in.close(), "CM input failed: "+path);}
	}
	private CovarianceModelParser(String path_, ByteFile in_){path=path_;in=in_;}

	private CovarianceModel parse(){
		assert(lineNumber==0):"Parser instances consume exactly one model";
		needLine();check(lp.termEquals("INFERNAL1/a", 0), "Only INFERNAL1/a is supported");
		final ObjectMap<String,String> headers=new ObjectMap<String,String>(String.class, String.class);
		while(true){needLine();if(lp.termEquals("CM", 0)){fields(1);break;}
			final String key=lp.parseString(0);check(lp.terms()>1, "Empty header "+key);
			check(knownHeader(key), "Unrecognized header "+key);
			lp.setBounds(1);final String value=new String(line, lp.a(), line.length-lp.a(), StandardCharsets.US_ASCII);
			if(key.equals("COM")){final String old=headers.get(key);headers.put(key, old==null ? value : old+"\n"+value);continue;}
			check(headers.put(key, value)==null, "Duplicate header "+key);
			validateHeader(key);
		}
		for(String key:new String[]{"NAME", "STATES", "NODES", "CLEN", "W", "ALPH", "NULL", "ELSELF", "EFP7GF"}){
			check(headers.get(key)!=null, "Missing required header "+key);}
		final int states=positiveHeader(headers, "STATES"), nodes=positiveHeader(headers, "NODES");
		final int clen=positiveHeader(headers, "CLEN"), window=positiveHeader(headers, "W");
		check(nodes<=states, "Each node requires at least one state");
		final byte[] nullBytes=headers.get("NULL").getBytes(StandardCharsets.US_ASCII);lp.set(nullBytes);
		final float[] nullModel=new float[4];for(int i=0; i<4; i++){nullModel[i]=prob(i, 0.25f);check(nullModel[i]>0, "NULL probabilities must be positive for log odds");}
		model=new CovarianceModel(headers.get("NAME"), headers.get("ACC"), clen, window, states, nodes, headers, nullModel);
		Arrays.fill(model.nodeState, -1);int nd=-1, nNodes=0;
		for(int v=0; v<states; v++){
			needLine();if(lp.termEquals("[", 0)){nd=parseNode(v, headers);nNodes++;needLine();}
			check(nd>=0, "A node must precede its states");parseState(v, nd);
		}
		needLine();fields(1);check(lp.termEquals("//", 0), "Missing closing CM //");
		check(nNodes==nodes, "Node count disagrees with NODES");
		// read_asc_1p1_cm(read_fp7=false) stops here. Validate the optional HMM envelope,
		// retaining no filter fields; this parser never claims to parse or use the HMM.
		if(nextLine()){
			check(lp.termStartsWith("HMMER3/", 0), "Expected optional filter HMM after CM");
			boolean closed=false;while(nextLine()){if(lp.termEquals("//", 0)){fields(1);closed=true;break;}}
			check(closed, "Unterminated optional filter HMM");check(!nextLine(), "Multiple models or trailing content unsupported");
		}
		validateTopology();normalize(model.nullModel);
		for(int v=0; v<states; v++){
			normalize(model.transition[v]);normalize(model.emission[v]);
			model.transitionScore[v]=scores(model.transition[v], false);
			model.emissionScore[v]=scores(model.emission[v], true);
		}
		consensusMap();return model;
	}

	private boolean knownHeader(String key){
		assert(key!=null):"Header dispatch must not conflate missing and unknown fields";
		for(String k:HEADERS){if(key.equals(k)){return true;}}return false;
	}
	private void validateHeader(String key){
		assert(lp.terms()>1):"Header values were checked before type validation";
		if(key.equals("DESC") || key.equals("DATE")){return;}
		if(key.equals("NULL")){fields(5);for(int i=1; i<5; i++){prob(i, 0.25f);}return;}
		if(key.equals("EFP7GF")){fields(3);number(1);number(2);return;}
		if(key.startsWith("ECM")){fields(7);for(int i=1; i<7; i++){number(i);}return;}
		fields(2);
		if(key.equals("NAME") || key.equals("ACC")){return;}
		if(key.equals("ALPH")){check(lp.termEquals("RNA", 1), "Only RNA alphabet supported");return;}
		if(key.equals("RF") || key.equals("CONS") || key.equals("MAP")){
			check(lp.termEquals("yes", 1) || lp.termEquals("no", 1), "Expected yes/no for "+key);return;}
		if(key.equals("STATES") || key.equals("NODES") || key.equals("CLEN") || key.equals("W") || key.equals("NSEQ")){
			check(integer(1)>0, "Positive integer required for "+key);return;}
		number(1);
	}
	private int positiveHeader(ObjectMap<String,String> h, String key){
		assert(h.get(key)!=null):"Required CM dimensions were checked before allocation";
		try{final int n=Integer.parseInt(h.get(key));check(n>0, "Positive "+key+" required");return n;}
		catch(NumberFormatException e){throw error("Invalid integer "+key);}
	}
	private int parseNode(int v, ObjectMap<String,String> headers){
		assert(v>=0 && v<model.states()):"Node mapping must point into the allocated state table";
		fields(10);check(lp.termEquals("]", 3), "Node closing bracket missing");
		final byte type=code(1, NODE_NAMES);final int n=integer(2);
		check(n>=0 && n<model.nodes() && model.nodeState[n]<0, "Invalid or duplicate node index "+n);
		model.nodeType[n]=type;model.nodeState[n]=v;
		final boolean left=type==MATP || type==MATL, right=type==MATP || type==MATR;
		model.mapLeft[n]=annotationMap(4, left && "yes".equals(headers.get("MAP")));
		model.mapRight[n]=annotationMap(5, right && "yes".equals(headers.get("MAP")));
		model.consLeft[n]=annotation(6, left && "yes".equals(headers.get("CONS")));
		model.consRight[n]=annotation(7, right && "yes".equals(headers.get("CONS")));
		model.rfLeft[n]=annotation(8, left && "yes".equals(headers.get("RF")));
		model.rfRight[n]=annotation(9, right && "yes".equals(headers.get("RF")));return n;
	}
	private int annotationMap(int field, boolean present){
		assert(field==4 || field==5):"CM node map occupies the two fields after the bracket";
		if(!present){check(lp.termEquals("-", field), "Unexpected node MAP annotation");return -1;}
		final int n=integer(field);check(n>0, "MAP positions are one-based");return n;
	}
	private byte annotation(int field, boolean present){
		assert(field>=6 && field<=9):"Node sequence/RF annotations use four single-character fields";
		check(lp.length(field)==1 && (present || lp.termEquals("-", field)), "Invalid node character annotation");return lp.parseByte(field, 0);
	}
	private void parseState(int v, int nd){
		assert(nd>=0):"Every state belongs to the preceding node record";
		check(lp.terms()>=10, "State requires ten structural fields");final byte t=code(0, STATE_NAMES);
		check(integer(1)==v, "Nonsequential state index; expected "+v);
		model.type[v]=t;model.node[v]=nd;model.parentLast[v]=integer(2);model.parentCount[v]=integer(3);
		model.childFirst[v]=integer(4);model.childCount[v]=integer(5);
		for(int i=0; i<4; i++){model.qdb[v][i]=integer(6+i);}
		final int nt=t==B ? 0 : model.childCount[v], ne=t==MP ? 16 : (t==ML || t==MR || t==IL || t==IR ? 4 : 0);
		check(nt>=0 && nt<=6, "Ordinary CM states have at most six children (Infernal MAXCONNECT)");fields(10+nt+ne);
		model.transition[v]=new float[nt];model.emission[v]=new float[ne];
		for(int i=0; i<nt; i++){model.transition[v][i]=prob(10+i, 1);}
		for(int i=0; i<ne; i++){final float q=ne==16 ? model.nullModel[i/4]*model.nullModel[i%4] : model.nullModel[i];model.emission[v][i]=prob(10+nt+i, q);}
	}
	private void validateTopology(){
		assert(model!=null):"Topology validation follows complete CM parsing";
		check(model.nodeType[0]==ROOT && model.nodeState[0]==0 && model.type[0]==S, "CM must begin at ROOT S");
		for(int n=0; n<model.nodes(); n++){check(model.nodeState[n]>=0, "Missing node "+n);}
		int clen=0;
		for(byte t:model.nodeType){clen+=(t==MATP ? 2 : (t==MATL || t==MATR ? 1 : 0));}
		check(clen==model.clen, "Node emissions disagree with declared CLEN");
		for(int v=0; v<model.states(); v++){
			final byte t=model.type[v];final int first=model.childFirst[v], count=model.childCount[v];
			check(model.parentCount[v]>=0 && model.parentLast[v]>=-1 && model.parentLast[v]<model.states(), "Invalid parent reference at state "+v);
			if(t==B){check(first>v && count>v && first<model.states() && count<model.states(), "B stores two descendant indices, not a child range");
				check(model.nodeType[model.node[first]]==BEGL && model.nodeType[model.node[count]]==BEGR, "B children must start left/right subtrees");}
			else if(t==E || t==EL){check(count==0, "Terminal state must not have children");}
			else{check(count>0 && first>=v && (long)first+count<=model.states(), "Invalid child interval at state "+v);
				check(first!=v || t==IL || t==IR, "Only consuming insertion states may self-loop");}
		}
	}
	private void consensusMap(){
		assert(model.nodes()>0):"Consensus traversal requires the validated root node";
		final IntList stack=new IntList();final boolean[] seen=new boolean[model.nodes()];stack.add(0);int cpos=0;
		while(stack.size>0){final int item=stack.pop(), n=item>>>1;check(n<model.nodes(), "Consensus traversal leaves node table");
			final byte t=model.nodeType[n];if((item&1)==1){model.consensusRight[n]=cpos+1;if(t==MATP || t==MATR){cpos++;}continue;}
			check(!seen[n], "Node topology revisits a subtree");seen[n]=true;
			if(t==MATP || t==MATL){cpos++;}model.consensusLeft[n]=cpos;stack.add((n<<1)|1);
			if(t==BIF){final int v=model.nodeState[n];check(model.type[v]==B, "BIF node must start at B state");
				stack.add(model.node[model.childCount[v]]<<1);stack.add(model.node[model.childFirst[v]]<<1);}
			else if(t!=END){stack.add((n+1)<<1);}
		}
		check(cpos==model.clen, "Consensus traversal must reproduce CLEN");for(boolean s:seen){check(s, "Unreachable node in CM topology");}
	}
	private void normalize(float[] p){
		assert(p!=null):"Every state has an explicit distribution, possibly empty";
		if(p.length==0){return;}float sum=0;for(float x:p){sum+=x;}
		check(Float.isFinite(sum) && sum>0, "Probability distribution has no finite positive mass");
		for(int i=0; i<p.length; i++){p[i]/=sum;}
	}
	private float[] scores(float[] p, boolean emit){
		assert(!emit || p.length==0 || p.length==4 || p.length==16):"Emission dimensions follow CM state type";
		final float[] s=new float[p.length];for(int i=0; i<s.length; i++){
			final float q=!emit ? 1 : (p.length==16 ? model.nullModel[i/4]*model.nullModel[i%4] : model.nullModel[i]);
			s[i]=p[i]==0 ? Float.NEGATIVE_INFINITY : (float)(Math.log(p[i]/q)/Math.log(2));}return s;
	}
	private float prob(int field, float nullFactor){
		assert(nullFactor>0):"Log-odds reconstruction requires positive null mass";
		if(lp.termEquals("*", field)){return 0;}final float p=(float)(Math.exp(number(field)/1.44269504)*nullFactor);
		check(Float.isFinite(p) && p>=0, "Invalid reconstructed probability");return p;
	}
	private double number(int field){
		assert(field>=0 && field<lp.terms()):"Numeric field must exist before parsing";
		// Header exponents need Java parsing; ordinary state scores stay on bytes.
		lp.setBounds(field);final byte[] data=lp.line();final int end=lp.b();int i=lp.a(), digits=0;
		if(i<end && (data[i]=='-' || data[i]=='+')){i++;}
		while(i<end && digit(data[i])){i++;digits++;}
		if(i<end && data[i]=='.'){i++;while(i<end && digit(data[i])){i++;digits++;}}
		check(digits>0, "Numeric field contains no digits");boolean exponent=false;
		if(i<end && (data[i]=='e' || data[i]=='E')){exponent=true;i++;if(i<end && (data[i]=='-' || data[i]=='+')){i++;}
			final int start=i;while(i<end && digit(data[i])){i++;}check(i>start, "Exponent contains no digits");}
		check(i==end, "Malformed numeric field "+field);
		try{final double x=exponent || digits>18 || data[lp.a()]=='+' ? Double.parseDouble(lp.parseString(field)) : lp.parseDouble(field);check(Double.isFinite(x), "Nonfinite number");return x;}
		catch(NumberFormatException e){throw error("Invalid numeric field "+field);}
	}
	private int integer(int field){
		assert(field>=0 && field<lp.terms()):"Integer field must exist before parsing";
		lp.setBounds(field);int a=lp.a();final int end=lp.b();final boolean negative=lp.line()[a]=='-';if(negative){a++;}
		check(a<end, "Integer field contains no digits");long x=0;
		for(; a<end; a++){final byte c=lp.line()[a];check(c>='0' && c<='9', "Invalid integer digit");x=x*10+c-'0';check(x<=2147483648L, "Integer overflow");}
		x=negative ? -x : x;check(x>=Integer.MIN_VALUE && x<=Integer.MAX_VALUE, "Integer overflow");return (int)x;
	}
	private static boolean digit(byte c){return c>='0' && c<='9';}
	private byte code(int field, String[] names){
		assert(names.length<128):"State/node codes are stored as nonnegative bytes";
		for(byte i=0; i<names.length; i++){if(lp.termEquals(names[i], field)){return i;}}throw error("Unknown state/node type: "+lp.parseString(field));
	}
	private boolean nextLine(){
		assert(in!=null):"CM parser must own an open ByteFile";
		while((line=in.nextLine())!=null){lineNumber++;int out=0;boolean space=false;
			for(int i=0; i<line.length; i++){final byte c=line[i];if(c=='#'){break;}
				if(c==' ' || c=='\t' || c=='\r'){space=out>0;}else{if(space){line[out++]=' ';space=false;}line[out++]=c;}}
			if(out==0){continue;}line=Arrays.copyOf(line, out);lp.set(line);return true;}
		return false;
	}
	private void needLine(){check(nextLine(), "Unexpected end of CM");}
	private void fields(int expected){check(lp.terms()==expected, "Expected "+expected+" fields, observed "+lp.terms());}
	private void check(boolean ok, String why){if(!ok){throw error(why);}}
	private IllegalArgumentException error(String why){return new IllegalArgumentException(path+":"+lineNumber+": "+why);}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	private final String path;
	private final ByteFile in;
	private final LineParser1 lp=new LineParser1(' ');
	private byte[] line;
	private int lineNumber;
	private CovarianceModel model;
	private static final String[] HEADERS={"NAME", "ACC", "DESC", "STATES", "NODES", "CLEN", "W", "ALPH", "RF", "CONS", "MAP", "DATE", "COM", "PBEGIN", "PEND", "WBETA", "QDBBETA1", "QDBBETA2", "N2OMEGA", "N3OMEGA", "ELSELF", "NSEQ", "EFFN", "CKSUM", "NULL", "GA", "TC", "NC", "EFP7GF", "ECMLC", "ECMGC", "ECMLI", "ECMGI"};
}
