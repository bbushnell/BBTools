package assemble;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import structures.ByteBuilder;

/** Optional visual attributes for Tadpole graphs; plain assembly DOT stays unchanged.
 * @author Fischl */
final class TadpoleDot {

	/** Opens a checked DOT output and writes its graph header. */
	static ByteStreamWriter open(final String path, final boolean pretty, final int k){
		assert(path!=null) : "A DOT destination is required before opening its writer.";
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TEXT, null, true,
				Tadpole.overwrite, Tadpole.append, false);
		final ByteStreamWriter writer=new ByteStreamWriter(ff);
		writer.start();
		writer.print("digraph G {\n");
		if(pretty){
			writer.print("\tgraph [overlap=false, rankdir=LR, nodesep=0.5, ranksep=1, pad=0.2, fontname=Helvetica, fontsize=11, labelloc=t, label=\"k="+
					k+"   |   depth=min,mean,max\"];\n");
			writer.print("\tnode [shape=box, fixedsize=false, style=\"rounded,filled\", margin=\"0.16,0.12\", fillcolor=\"#e8f0fa\", color=\"#677d96\", fontname=Helvetica, fontsize=11];\n");
			writer.print("\tedge [fontname=Helvetica, fontsize=10, color=\"#46566b\", arrowsize=0.65];\n");
		}
		return writer;
	}

	/** Appends one node; external names are escaped, never used as structural IDs. */
	static void node(final Contig c, final int offset, final int k, final boolean pretty,
			final boolean graphOnly, final double logLengthPerInch, final ByteBuilder bb){
		assert(c!=null && c.bases!=null) : "Graph nodes must retain their source sequence.";
		bb.tab().append(c.id).append(" [label=\"");
		if(pretty || c.name!=null){
			if(c.name!=null){escape(c.name, bb);}
			else{bb.append("contig_").append(c.id+offset);}
			bb.append("\\n");
		}
		if(!pretty){bb.append("id=").append(c.id).append("\\n");}
		bb.append("len=").append(c.length());
		if(pretty){
			bb.append("\\ndepth=").append(c.minCov).append(',').append(c.coverage, 2).append(',').append(c.maxCov);
		}else{bb.append("\\ncov=").append(c.coverage, 1);}
		if(graphOnly && !pretty){bb.append("\\nk=").append(k);}
		if(!graphOnly){
			bb.append("\\nleft=").append(Tadpole.codeStrings[c.leftCode]);
			bb.append("\\nright=").append(Tadpole.codeStrings[c.rightCode]);
			if(c.graphClass>=0){bb.append("\\nclass="); c.appendGraphClass(bb);}
		}
		if(c.graphLeftLimited || c.graphRightLimited){
			bb.append("\\nsearch_limited=");
			if(c.graphLeftLimited){bb.append('L');}
			if(c.graphRightLimited){bb.append('R');}
		}
		bb.append('"');
		if(pretty){
			//Logarithmic widths keep short contigs visible. Graphviz may enlarge a
			//node further to fit its original name; labels must remain inside it.
			bb.append(", width=").append(nodeWidth(c.length(), logLengthPerInch), 4);
			bb.append(", height=0.75");
		}
		bb.append("]\n");
	}

	/** Normalizes log(max(16,length-128)) so the longest contig requests four inches. */
	static double logLengthPerInch(final ArrayList<Contig> contigs){
		assert(contigs!=null) : "DOT scaling requires the actual graph's contig lengths.";
		int max=1;
		for(Contig c : contigs){max=Math.max(max, c.length());}
		return Math.log(Math.max(16, max-128))/4.0;
	}

	/** Minimum requested width; label fitting can expand it, never clip the name. */
	static double nodeWidth(final int length, final double logLengthPerInch){
		assert(length>=0 && logLengthPerInch>0) : "DOT widths require a nonnegative length and positive log scale.";
		return Math.max(1.2, Math.log(Math.max(16, length-128))/logLengthPerInch);
	}

	/** Appends optional depth styling without changing the legacy edge representation. */
	static void edge(final Edge e, final int k, final boolean pretty, final ByteBuilder bb){
		assert(e!=null) : "DOT edge serialization requires a connection.";
		if(!pretty && e.pathMinDepth<0){e.toDot(bb); return;}
		bb.append(e.origin).append(" -> ").append(e.destination).append(" [label=\"");
		if(pretty){
			bb.append(e.sourceRight() ? 'R' : 'L').append("->").append(e.destRight() ? 'R' : 'L');
			bb.append(" steps=").append(e.length);
			if(e.pathMinDepth>=0){
				assert(e.pathMaxDepth>=e.pathMinDepth) : "A displayed path depth range must come from the same measured route.";
				bb.append("\\ndepth=").append(e.pathMinDepth).append(',').append(e.pathMeanDepth, 2).append(',').append(e.pathMaxDepth);
			}else{bb.append("\\nstartdepth=").append(e.depth);}
		}else{
			bb.append(e.sourceRight() ? "RIGHT" : "LEFT");
			bb.append("\\nlen=").append(e.length).append("\\norient=").append(e.orientation);
			bb.append("\\nk=").append(k).append("\\nstartdepth=").append(e.depth);
		}
		if(!pretty && e.pathMinDepth>=0){
			bb.append("\\nmin=").append(e.pathMinDepth);
			bb.append("\\nmean=").append(e.pathMeanDepth, 2);
		}
		bb.append('"');
		if(pretty){
			final int depth=(e.pathMinDepth>=0 ? e.pathMinDepth : e.depth);
			bb.append(", penwidth=").append(penWidth(depth), 3);
			bb.append(", tailport=").append(e.sourceRight() ? 'e' : 'w');
			bb.append(", headport=").append(e.destRight() ? 'e' : 'w');
		}
		bb.append("]\n");
	}

	/** Compresses the depth range visually; the exact depth stays in the edge label. */
	static double penWidth(final int depth){
		assert(depth>=0) : "Read kmer support cannot be negative in a graph label: "+depth;
		return 0.65+Math.min(0.65, 0.08*Math.log1p(Math.max(0, depth))/Math.log(2));
	}

	/** Escapes DOT quoted text, including literal backslashes and line breaks. */
	static void escape(final String text, final ByteBuilder bb){
		assert(text!=null && bb!=null) : "DOT escaping needs source text and an output buffer.";
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);
			if(c=='\\' || c=='"'){bb.append('\\').append(c);}
			else if(c=='\n'){bb.append("\\n");}
			else if(c=='\r'){bb.append("\\r");}
			else if(c=='\t'){bb.append("\\t");}
			else{bb.append(c);}
		}
	}
}
