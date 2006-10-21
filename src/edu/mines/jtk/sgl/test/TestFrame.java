/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package edu.mines.jtk.sgl.test;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;

import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.sgl.*;

/**
 * A frame for testing the package {@link edu.mines.jtk.sgl}.
 * @author Dave Hale
 * @version 2006.06.28
 */
public class TestFrame extends JFrame {
  private static final long serialVersionUID = 1L;

  public TestFrame(World world) {
    OrbitView view = (world!=null)?new OrbitView(world):new OrbitView();
    view.setAxesOrientation(View.AxesOrientation.XRIGHT_YOUT_ZDOWN);
    ViewCanvas canvas = new ViewCanvas(view);
    canvas.setView(view);

    ModeManager mm = new ModeManager();
    mm.add(canvas);
    OrbitViewMode ovm = new OrbitViewMode(mm);
    SelectDragMode sdm = new SelectDragMode(mm);

    JPopupMenu.setDefaultLightWeightPopupEnabled(false);
    ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);

    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic('F');
    Action exitAction = new AbstractAction("Exit") {
      private static final long serialVersionUID = 1L;
      public void actionPerformed(ActionEvent event) {
        System.exit(0);
      }
    };
    JMenuItem exitItem = fileMenu.add(exitAction);
    exitItem.setMnemonic('x');

    JMenu modeMenu = new JMenu("Mode");
    modeMenu.setMnemonic('M');
    JMenuItem ovmItem = new ModeMenuItem(ovm);
    modeMenu.add(ovmItem);
    JMenuItem sdmItem = new ModeMenuItem(sdm);
    modeMenu.add(sdmItem);

    JMenuBar menuBar = new JMenuBar();
    menuBar.add(fileMenu);
    menuBar.add(modeMenu);

    JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
    toolBar.setRollover(true);
    JToggleButton ovmButton = new ModeToggleButton(ovm);
    toolBar.add(ovmButton);
    JToggleButton sdmButton = new ModeToggleButton(sdm);
    toolBar.add(sdmButton);

    ovm.setActive(true);

    this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    this.setSize(new Dimension(SIZE,SIZE));
    this.add(canvas,BorderLayout.CENTER);
    this.add(toolBar,BorderLayout.WEST);
    this.setJMenuBar(menuBar);

    _view = view;
  }

  public OrbitView getOrbitView() {
    return _view;
  }

  public static void main(String[] args) {
    TestFrame frame = new TestFrame(null);
    frame.setVisible(true);
  }

  private static final int SIZE = 600;
  private OrbitView _view;
}
