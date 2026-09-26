package gff;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Tools;

/**
 * Reusable, allocation-light checks for GFF core fields.
 * Checks nine nonempty tab-separated fields, unescaped control characters,
 * positive ordered coordinates, decimal score syntax, strand, and phase.
 * Coordinates are compared as decimal text, without an artificial int/long limit.
 * Attribute syntax, directives, ontology membership, feature relationships, and
 * embedded FASTA are outside this check's scope. Passing is not full GFF3 validation.
 * Rules: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
 * @author Shinobu
 * @date Sep 25, 2026
 */
public final class GffValidator {

	private GffValidator(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks one feature line using a fresh parser. For repeated calls, reuse the
	 * scratch-parser overload to avoid allocating a parser per record.
	 * @param line Feature bytes without a line terminator; comments are not features
	 * @return Null on success, otherwise the first core-field diagnostic
	 */
	public static String validateFields(final byte[] line){
		return validateFields(line, new LineParser1('\t'));
	}

	/**
	 * Checks one feature line without constructing field Strings or GffLine objects.
	 * The parser is overwritten and must belong exclusively to this invocation or
	 * worker. The input is not modified; separate parsers allow concurrent calls.
	 * Malformed fields return diagnostics even when assertions are disabled.
	 * @param line Feature bytes without CR/LF; null is reported as invalid
	 * @param parser Caller-owned tab parser, reused across records
	 * @return Null if the documented core checks pass, otherwise a diagnostic
	 * @throws IllegalArgumentException If scratch is null or uses another delimiter
	 */
	public static String validateFields(final byte[] line, final LineParser1 parser){
		if(parser==null || parser.delimiter!='\t'){
			throw new IllegalArgumentException("GFF field validation requires a caller-owned tab parser.");
		}
		if(line==null || line.length==0){return "Missing feature line";}
		parser.set(line);
		if(parser.terms()!=9){return "Expected 9 tab-separated fields; found "+parser.terms();}
		for(int field=0; field<9; field++){
			if(parser.length(field)==0){return "Empty field "+(field+1)+"; use '.' for an undefined value";}
		}
		for(byte b : line){
			if((b>=0 && b<32 && b!='\t') || b==127){return "Unescaped control character in feature line";}
		}
		parser.setBounds(3);
		final int startA=parser.a(), startB=parser.b();
		parser.setBounds(4);
		final int endA=parser.a(), endB=parser.b();
		if(!positiveInteger(line, startA, startB)){return "Field 4 (start) must be a positive decimal integer";}
		if(!positiveInteger(line, endA, endB)){return "Field 5 (end) must be a positive decimal integer";}
		if(compareCoordinates(line, startA, startB, endA, endB)>0){return "Field 4 (start) exceeds field 5 (end)";}
		parser.setBounds(5);
		if(!parser.currentTermEquals((byte)'.') && !decimal(line, parser.a(), parser.b())){
			return "Field 6 (score) must be a decimal number, optionally with an exponent, or '.'";
		}
		if(!(parser.termEquals('+', 6) || parser.termEquals('-', 6) ||
			parser.termEquals('.', 6) || parser.termEquals('?', 6))){return "Field 7 (strand) must be '+', '-', '.', or '?'";}
		if(!(parser.termEquals('.', 7) || parser.termEquals('0', 7) ||
			parser.termEquals('1', 7) || parser.termEquals('2', 7))){return "Field 8 (phase) must be '.', 0, 1, or 2";}
		if((parser.termEquals("CDS", 2) || parser.termEquals("SO:0000316", 2)) && parser.termEquals('.', 7)){
			return "CDS features require phase 0, 1, or 2";
		}
		return null;
	}

	/**
	 * Scans a file to EOF with private scratch and a single-threaded byte reader.
	 * Blank lines and comments are counted but not validated. An exact ##FASTA
	 * directive or an implicit '>' header starts an unvalidated tail (the GFF3
	 * Artemis compatibility rule), which is still consumed to detect I/O
	 * failures. At most one feature diagnostic is retained; all features are counted.
	 * @param ff Input format, including compression and subprocess policy
	 * @return Core-check counts, first diagnostic, and latched read/close error
	 */
	public static Result validateFile(final FileFormat ff){
		if(ff==null || !ff.read()){throw new IllegalArgumentException("GFF validation requires an input FileFormat.");}
		final ByteFile bf=new ByteFile1(ff);
		final LineParser1 parser=new LineParser1('\t');
		long lines=0, features=0, invalid=0, comments=0, blanks=0, fastaLines=0, firstErrorLine=0;
		String firstError=null;
		boolean inFasta=false, readError;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lines++;
				if(inFasta){fastaLines++; continue;}
				if(line.length==0){blanks++; continue;}
				if(line[0]=='>'){inFasta=true; fastaLines++; continue;}
				if(line[0]=='#'){
					comments++;
					if(line.length==7 && Tools.startsWith(line, "##FASTA")){inFasta=true;}
					continue;
				}
				features++;
				final String error=validateFields(line, parser);
				if(error!=null){
					invalid++;
					if(firstError==null){firstError=error; firstErrorLine=lines;}
				}
			}
		}finally{
			readError=bf.close();
		}
		assert(lines==features+comments+blanks+fastaLines) : "GFF scan must classify each consumed line exactly once";
		return new Result(lines, features, invalid, comments, blanks, fastaLines, firstErrorLine, firstError, readError);
	}

	/** Tests positive decimal syntax without overflow; leading zeroes are allowed. */
	private static boolean positiveInteger(final byte[] line, final int a, final int b){
		assert(a>=0 && a<b && b<=line.length) : "Nonempty field bounds are established by validateFields";
		boolean nonzero=false;
		for(int i=a; i<b; i++){
			final byte c=line[i];
			if(c<'0' || c>'9'){return false;}
			nonzero|=c!='0';
		}
		return nonzero;
	}

	/** Compares already validated positive integers after discarding leading zeroes. */
	private static int compareCoordinates(final byte[] line, int a, final int b, int c, final int d){
		assert(a<b && c<d) : "Coordinate comparison requires nonempty positiveInteger-validated fields";
		while(line[a]=='0'){a++;}
		while(line[c]=='0'){c++;}
		if(b-a!=d-c){return b-a<d-c ? -1 : 1;}
		while(a<b){
			if(line[a]!=line[c]){return line[a]-line[c];}
			a++; c++;
		}
		return 0;
	}

	/** Decimal lexical check, including signs/exponents; excludes NaN, infinity and hex. */
	private static boolean decimal(final byte[] line, int a, final int b){
		assert(a>=0 && a<b && b<=line.length) : "Score bounds must describe the nonempty sixth GFF field";
		if(line[a]=='+' || line[a]=='-'){a++;}
		int digits=0;
		while(a<b && line[a]>='0' && line[a]<='9'){a++; digits++;}
		if(a<b && line[a]=='.'){
			a++;
			while(a<b && line[a]>='0' && line[a]<='9'){a++; digits++;}
		}
		if(digits==0){return false;}
		if(a<b && (line[a]=='e' || line[a]=='E')){
			a++;
			if(a<b && (line[a]=='+' || line[a]=='-')){a++;}
			final int exponentStart=a;
			while(a<b && line[a]>='0' && line[a]<='9'){a++;}
			if(a==exponentStart){return false;}
		}
		return a==b;
	}

	/** Immutable result. A passing core-field scan does not imply complete GFF3 validity. */
	public static final class Result {
		private Result(final long lines_, final long features_, final long invalid_, final long comments_,
			final long blanks_, final long fastaLines_, final long firstErrorLine_, final String firstError_, final boolean readError_){
			lines=lines_; features=features_; invalid=invalid_; comments=comments_; blanks=blanks_;
			fastaLines=fastaLines_; firstErrorLine=firstErrorLine_; firstError=firstError_; readError=readError_;
		}
		/** True only if all visited feature lines passed the core checks and reading succeeded. */
		public boolean passedCoreChecks(){return invalid==0 && !readError;}
		/** Total consumed physical lines and feature/comment/blank/unvalidated-tail counts. */
		public final long lines, features, invalid, comments, blanks, fastaLines;
		/** One-based first invalid feature line, or zero when no feature error was found. */
		public final long firstErrorLine;
		/** First core-field diagnostic, or null; does not include file contents. */
		public final String firstError;
		/** Latched read/decompression/close failure reported by the byte reader. */
		public final boolean readError;
	}
}
