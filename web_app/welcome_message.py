from dash import dcc, html
from styling import text_style

welcome_message_content = html.Div([
    dcc.Markdown(
        """
        ## **Introduction to StreptoCAD**

        
        - StreptoCAD is an open-source software toolbox designed to help you build biology easier.
        It supports genome engineering workflows in Streptomyces, providing a user-friendly interface to design your experiments.
        We hope you find StreptoCAD useful in your research endeavors. If you have any questions or feedback, please don't hesitate to reach out to our team.


        ### **Getting Started**

        
        - 👈 To get started, please select a workflow from the tabs on the left.

        - If you need help in any way, you are welcome to contact us. You can find contact details in the bottom the page.

        ### **Happy bioengineering!**

        """,
        style={**text_style, "marginBottom": "200px", "lineHeight": "2", "marginTop": "30px"}
    )])
