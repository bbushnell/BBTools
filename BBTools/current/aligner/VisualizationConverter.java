package aligner;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import javax.imageio.ImageIO;

/**
 * Converts text-based alignment visualizations to color bitmap images.
 * Uses a color gradient to represent alignment scores with:
 * - Green (a-z): Low scores
 * - Yellow (0-9): Medium scores
 * - Red (A-Z): High scores
 * - Black: Unexplored regions
 * - Gray: Explored but pruned regions
 * - White: Optimal alignment path
 * @author Brian Bushnell
 * @author Isla (Highly-customized Claude instance)
 * @date April 2025
 */
public class VisualizationConverter {
	// Pixel scale factor - makes each character larger in the output image
	private static final int DEFAULT_SCALE = 2;

	public static void main(String[] args) {
		if(args.length==1 && args[0].endsWith(".txt")) {
			args=new String[] {args[0], args[0].replace(".txt", ".png")};
		}
		if(args.length>2) {colorScheme=Integer.parseInt(args[2]);}
		if (args.length < 2) {
			System.err.println("Usage: java aligner.VisualizationConverter <input_text_file> <output_image>");
			System.exit(1);
		}

		String textFile = args[0];
		String imageFile = args[1];

		try {
			convertToBitmap(textFile, imageFile);
			System.out.println("Successfully converted " + textFile + " to " + imageFile);
		} catch (IOException e) {
			System.err.println("Error during conversion: " + e.getMessage());
			e.printStackTrace();
		}
	}

	private static void convertToBitmap(String textFile, String imageFile) throws IOException {
		// Read text file to determine dimensions
		int width = 0;
		int height = 0;

		try (BufferedReader reader = new BufferedReader(new FileReader(textFile))) {
			String line;
			while ((line = reader.readLine()) != null) {
				width = Math.max(width, line.length());
				height++;
			}
		}

		if (width == 0 || height == 0) {
			throw new IOException("Input file is empty or could not be read");
		}
		int scale=(width>=2000 ? 1 : DEFAULT_SCALE);

		// Create the image with scaling
		BufferedImage image = new BufferedImage(width * scale, height * scale, BufferedImage.TYPE_INT_RGB);

		// Fill the image with data from the text file
		try (BufferedReader reader = new BufferedReader(new FileReader(textFile))) {
			String line;
			int y = 0;

			while ((line = reader.readLine()) != null && y < height) {
				for (int x = 0; x < line.length(); x++) {
					char c = line.charAt(x);
					Color color = (colorScheme==1 ? getColorForChar1(c) : getColorForChar2(c));

					// Fill the scaled pixel block
					for (int sy = 0; sy < scale; sy++) {
						for (int sx = 0; sx < scale; sx++) {
							image.setRGB(x * scale + sx, y * scale + sy, color.getRGB());
						}
					}
				}
				y++;
			}
		}

		// Determine output format from file extension
		String format = "png"; // Default
		int dotIndex = imageFile.lastIndexOf('.');
		if (dotIndex > 0) {
			String extension = imageFile.substring(dotIndex + 1).toLowerCase();
			if (extension.equals("jpg") || extension.equals("jpeg") || 
					extension.equals("png") || extension.equals("bmp")) {
				format = extension.equals("jpeg") ? "jpg" : extension;
			}
		}

		// Write the image to file
		File outputFile = new File(imageFile);
		ImageIO.write(image, format, outputFile);
	}


	private static Color getColorForChar(char level0, int scheme) {
		return scheme==1 ? getColorForChar1(level0) : getColorForChar2(level0);
	}


	private static Color getColorForChar1(char level0) {
		switch (level0) {
			case ' ': return Color.BLACK;            // Unexplored
			case '.': return new Color(80, 80, 80);  // Explored but pruned
			case '*': return Color.WHITE;            // Optimal path
		}
		
		final int level=Visualizer.symbolMap[level0];
		
//		if (level <= 15) {// Green range (16 steps)
//			float position = level / 15f;
//			return interpolateColor(new Color(0, 96, 0), new Color(0, 255, 0), position);
//		} else if (level <= 25) {// Green - Yellow range (10 steps)
//			float position = level / 9f;
//			return interpolateColor(new Color(0, 255, 0), new Color(240, 255, 0), position);
//		} else if (level <= 35) {// Yellow-orange range (10 steps)
//			float position = level / 9f;
//			return interpolateColor(new Color(255, 255, 0), new Color(255, 180, 0), position);
//		} else if (level <= 55) {// Orange-red (20 steps)
//			float position = level / 19f;
//			return interpolateColor(new Color(255, 170, 0), new Color(255, 16, 0), position);
//		} else if (level <= 62) {// Red-purple (6 steps) 
//			float position = level / 5f;
//			return interpolateColor(new Color(255, 0, 16), new Color(96, 32, 255), position);
//		}
		
		// Enhanced contrast for score characters (a-z, 0-9, A-Z)
		// Letter order: abcde fghij klmno pqrst uvwxy z
		if (level0 >= 'a' && level0 <= 'p') {
			// Green range (16 steps)
			float position = (float)(level0 - 'a') / 15f;
			return interpolateColor(new Color(0, 96, 0), new Color(0, 255, 0), position);
		} else if (level0 >= 'q' && level0 <= 'z') {
			// Green - Yellow range (10 steps)
			float position = (float)(level0 - 'q') / 9f;
			return interpolateColor(new Color(0, 255, 0), new Color(240, 255, 0), position);
		} else if (level0 >= '0' && level0 <= '9') {
			// Yellow-orange range (10 steps)
			float position = (float)(level0 - '0') / 9f;
			return interpolateColor(new Color(255, 255, 0), new Color(255, 180, 0), position);
		} else if (level0 >= 'A' && level0 <= 'T') {//20 steps
			// Orange-red (20 steps)
			float position = (float)(level0 - 'A') / 19f;
			return interpolateColor(new Color(255, 170, 0), new Color(255, 16, 0), position);
		} else if (level0 > 'T' && level0 <= 'Z') {
			// Red-purple (6 steps) 
			float position = (float)(level0 - 'T') / 5f;
			return interpolateColor(new Color(255, 0, 16), new Color(96, 32, 255), position);
		}

		// Default for any unrecognized character
		return Color.GRAY;
	}
	
	private static Color getColorForChar2(char level0) {
	    switch (level0) {
	        case ' ': return Color.BLACK;            // Unexplored
	        case '.': return new Color(80, 80, 80);  // Explored but pruned
	        case '*': return Color.WHITE;            // Optimal path
	    }
	    
	    final int level=Visualizer.symbolMap[level0];
	    
	    // Enhanced contrast for score characters (a-z, 0-9, A-Z)
	    if (level0 >= 'a' && level0 <= 'p') {
	        // Blue range (16 steps)
	        float position = (float)(level0 - 'a') / 15f;
	        return interpolateColor(new Color(0, 0, 128), new Color(0, 100, 255), position);
	    } else if (level0 >= 'q' && level0 <= 'z') {
	        // Blue - Cyan range (10 steps)
	        float position = (float)(level0 - 'q') / 9f;
	        return interpolateColor(new Color(0, 100, 255), new Color(0, 210, 210), position);
	    } else if (level0 >= '0' && level0 <= '9') {
	        // Cyan-Yellow range (10 steps)
	        float position = (float)(level0 - '0') / 9f;
	        return interpolateColor(new Color(0, 210, 210), new Color(240, 240, 0), position);
	    } else if (level0 >= 'A' && level0 <= 'T') {
	        // Yellow-Purple (20 steps)
	        float position = (float)(level0 - 'A') / 19f;
	        return interpolateColor(new Color(240, 240, 0), new Color(170, 0, 255), position);
	    } else if (level0 > 'T' && level0 <= 'Z') {
	        // Purple-Magenta (6 steps) 
	        float position = (float)(level0 - 'T') / 5f;
	        return interpolateColor(new Color(170, 0, 255), new Color(255, 0, 170), position);
	    }

	    // Default for any unrecognized character
	    return Color.GRAY;
	}

	private static Color interpolateColor(Color c1, Color c2, float position) {
		int red = (int)(c1.getRed() + position * (c2.getRed() - c1.getRed()));
		int green = (int)(c1.getGreen() + position * (c2.getGreen() - c1.getGreen()));
		int blue = (int)(c1.getBlue() + position * (c2.getBlue() - c1.getBlue()));
		return new Color(
				Math.max(0, Math.min(255, red)), 
				Math.max(0, Math.min(255, green)), 
				Math.max(0, Math.min(255, blue))
				);
	}
	
	public static int colorScheme=1;
	
}