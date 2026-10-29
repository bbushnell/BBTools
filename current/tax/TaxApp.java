package tax;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

/** Unfinished Swing demonstration containing a button and an input dialog.
 * It performs no taxonomy lookup and ignores the dialog result. This is not an
 * entry point for the working command-line taxonomy tools.
 * @author Brian Bushnell */
public class TaxApp{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Displays the placeholder window; clicking its button opens an input dialog.
	 * Window creation runs on the calling thread. No explicit close operation is
	 * configured, and the entered name is not retained or processed.
	 * @param args Unused command-line arguments */
	public static void main(final String[] args){
		final JFrame parent=new JFrame();
		JButton button=new JButton();
		button.setText("Button text");
		parent.add(button);
		parent.pack();
		parent.setVisible(true);
		button.addActionListener(new java.awt.event.ActionListener(){
			/** Opens the demonstration prompt and discards its result.
			 * @param evt Button action event, unused */
			@Override
			public void actionPerformed(java.awt.event.ActionEvent evt){
				String name=JOptionPane.showInputDialog(parent, "Prompt", null);
			}
		});
	}

	/** Creates an otherwise stateless demonstration holder; does not open a window. */
	public TaxApp(){}
}
