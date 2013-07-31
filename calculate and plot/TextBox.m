function [ ] = TextBox( text, position, mFigure )
%Create text box to figure
% text in string
% postion in [x y length height]
% mFigure is figure itself

% Create the figure
%mFigure = figure()

% Create a uicontrol of type "text"
mTextBox = uicontrol('style','text');
set(mTextBox,'String',text)

% Something that I find useful is to set the Position Units to Characters,the default is pixels
set(mTextBox,'Units','characters')
% This means a Text Box with 3 lines of text will have a height of 3

% To move the the Text Box around you can set and get the position of TextBox itself
mTextBoxPosition = get(mTextBox,'Position');
% The array mTextBoxPosition has four elements
% [x y length height]

%set textbox position
set(mTextBox,'Position',mTextBoxPosition+position)

% % Get the Color of the figure window
% colorOfFigureWindow = get(mFigure,'Color');
% 
% %Set the BackgroundColor of the text box
% set(mTextBox,'BackgroundColor',colorOfFigureWindow)

end

